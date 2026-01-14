"""Memory-efficient three-stage pipeline implementation.

This version minimizes memory usage by:
1. Using boolean masks instead of DataFrame copies
2. Pre-extracting numpy arrays for filter columns
3. Using numpy operations instead of pandas groupby during trials
4. Processing variant types sequentially with explicit memory cleanup
"""

import pandas as pd
import numpy as np
import optuna
import statsmodels.formula.api as smf
from typing import Dict, List, Any, Optional, Tuple
from joblib import Memory
import gc
import logging
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig
from .data import VariantDataset
from .plotting import plot_optimization_results, save_parameters


logger = logging.getLogger(__name__)


@dataclass
class OptimisationResult:
    """Container for optimisation results."""
    best_params: Dict[str, Dict[str, Any]]
    best_scores: Dict[str, float]
    warmup_params: Optional[Dict[str, Dict[str, Any]]]
    n_individuals_removed: int
    study: Optional[optuna.Study]
    
    @property
    def summary(self) -> str:
        """Human-readable summary."""
        lines = ["Optimisation Results:"]
        for var_type, params in self.best_params.items():
            lines.append(f"\n{var_type}:")
            score = self.best_scores.get(var_type, 0)
            score_str = f"{score:.4e}" if score < 0.01 else f"{score:.4f}"
            lines.append(f"  Best score: {score_str}")
            for param, value in params.items():
                if isinstance(value, float):
                    lines.append(f"  {param}: {value:.4f}")
                else:
                    lines.append(f"  {param}: {value}")
        if self.n_individuals_removed > 0:
            lines.append(f"\nOutliers removed: {self.n_individuals_removed}")
        return '\n'.join(lines)


@dataclass
class OptimisationArrays:
    """Pre-extracted numpy arrays for memory-efficient optimization.
    
    Instead of passing full DataFrames to the objective function,
    we pre-extract only the columns needed as numpy arrays.
    """
    # Sample IDs encoded as integers for fast groupby
    sample_ids: np.ndarray  # int array of sample indices
    sample_id_to_idx: Dict[str, int]  # mapping from sample name to index
    idx_to_sample_id: Dict[int, str]  # reverse mapping
    n_samples: int
    
    # Filter columns as numpy arrays (column_name -> array)
    filter_arrays: Dict[str, np.ndarray]
    
    # Parental ages per sample (indexed by sample_idx)
    paternal_ages: np.ndarray  # shape (n_samples,)
    maternal_ages: np.ndarray  # shape (n_samples,)
    valid_age_mask: np.ndarray  # bool array for samples with valid ages
    
    # Original number of variants
    n_variants: int


def extract_optimisation_arrays(
    df: pd.DataFrame,
    filter_columns: List[str],
    sample_col: str = 'SAMPLE'
) -> OptimisationArrays:
    """Extract numpy arrays from DataFrame for memory-efficient optimization.
    
    This function extracts only the data needed for the objective function
    as compact numpy arrays, allowing the original DataFrame to be freed.
    """
    # Create sample ID mapping
    unique_samples = df[sample_col].unique()
    sample_id_to_idx = {s: i for i, s in enumerate(unique_samples)}
    idx_to_sample_id = {i: s for s, i in sample_id_to_idx.items()}
    n_samples = len(unique_samples)
    
    # Encode sample IDs as integers
    sample_ids = df[sample_col].map(sample_id_to_idx).values.astype(np.int32)
    
    # Extract filter columns
    filter_arrays = {}
    for col in filter_columns:
        if col in df.columns:
            # Convert to float32 to save memory
            filter_arrays[col] = pd.to_numeric(df[col], errors='coerce').values.astype(np.float32)
    
    # Extract parental ages per sample
    sample_ages = df[[sample_col, 'paternal_age', 'maternal_age']].drop_duplicates(subset=[sample_col])
    sample_ages = sample_ages.set_index(sample_col)
    
    paternal_ages = np.full(n_samples, np.nan, dtype=np.float32)
    maternal_ages = np.full(n_samples, np.nan, dtype=np.float32)
    
    for sample, idx in sample_id_to_idx.items():
        if sample in sample_ages.index:
            paternal_ages[idx] = sample_ages.loc[sample, 'paternal_age']
            maternal_ages[idx] = sample_ages.loc[sample, 'maternal_age']
    
    valid_age_mask = ~(np.isnan(paternal_ages) | np.isnan(maternal_ages))
    
    return OptimisationArrays(
        sample_ids=sample_ids,
        sample_id_to_idx=sample_id_to_idx,
        idx_to_sample_id=idx_to_sample_id,
        n_samples=n_samples,
        filter_arrays=filter_arrays,
        paternal_ages=paternal_ages,
        maternal_ages=maternal_ages,
        valid_age_mask=valid_age_mask,
        n_variants=len(df)
    )


def apply_filters_mask(
    arrays: OptimisationArrays,
    params: Dict[str, Any]
) -> np.ndarray:
    """Apply filter parameters and return boolean mask (no copy created).
    
    Returns a boolean array indicating which variants pass all filters.
    """
    mask = np.ones(arrays.n_variants, dtype=bool)
    
    for param, value in params.items():
        if param.startswith('min_'):
            col = param[4:]
            if col in arrays.filter_arrays:
                col_mask = arrays.filter_arrays[col] >= value
                # Handle NaN values - exclude them
                col_mask &= ~np.isnan(arrays.filter_arrays[col])
                mask &= col_mask
        elif param.startswith('max_'):
            col = param[4:]
            if col in arrays.filter_arrays:
                col_mask = arrays.filter_arrays[col] <= value
                col_mask &= ~np.isnan(arrays.filter_arrays[col])
                mask &= col_mask
    
    return mask


def count_per_sample_fast(
    sample_ids: np.ndarray,
    mask: np.ndarray,
    n_samples: int
) -> np.ndarray:
    """Fast per-sample counting using numpy bincount.
    
    Much faster than pandas groupby for this use case.
    """
    # Only count samples that pass the filter
    filtered_samples = sample_ids[mask]
    counts = np.bincount(filtered_samples, minlength=n_samples)
    return counts


class MemoryEfficientPipeline:
    """Memory-efficient three-stage optimisation pipeline."""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        
        # Set up caching
        if config.use_cache:
            cache_dir = config.cache_dir or './cache'
            self.memory = Memory(cache_dir, verbose=0)
        
        # Set seeds for reproducibility
        if config.deterministic:
            self._set_all_seeds(config.seed)
    
    def _set_all_seeds(self, seed: int):
        """Set all random seeds for deterministic results."""
        import random
        import os
        random.seed(seed)
        np.random.seed(seed)
        os.environ['PYTHONHASHSEED'] = str(seed)
        optuna.logging.set_verbosity(optuna.logging.INFO)
    
    def run(
        self,
        data: VariantDataset,
        reference: VariantDataset,
        output_dir: Optional[Path] = None,
        generate_plots: bool = True
    ) -> OptimisationResult:
        """Run complete three-stage optimisation pipeline."""
        logger.info("="*60)
        logger.info("Starting MEMORY-EFFICIENT three-stage optimisation pipeline")
        logger.info("="*60)
        
        # Log initial memory usage
        self._log_memory_usage("Initial")

        # Calculate regression targets from reference
        logger.info("Calculating regression targets from reference data...")
        targets_by_type = self._calculate_targets(reference)
        logger.info(f"Targets calculated for {len(targets_by_type)} variant types")

        warmup_params = {}
        n_removed = 0
        samples_to_keep = None

        # Stage 1: Warmup - process each variant type sequentially
        logger.info("")
        logger.info("="*60)
        logger.info("STAGE 1: WARMUP OPTIMISATION")
        logger.info("="*60)
        if self.config.stage1.enabled:
            logger.info(f"Running warmup with {self.config.stage1.n_trials} trials per variant type...")
            warmup_params = self._run_warmup_sequential(data, targets_by_type)
            logger.info(f"✓ Warmup complete for {len(warmup_params)} variant types")
            self._log_memory_usage("After warmup")
        else:
            logger.info("Stage 1 disabled in configuration")

        # Stage 2: Outlier removal
        logger.info("")
        logger.info("="*60)
        logger.info("STAGE 2: OUTLIER REMOVAL")
        logger.info("="*60)
        if self.config.stage2.enabled and warmup_params:
            logger.info(f"Removing outliers with DNM count range: {self.config.stage2.min_dnm_count}-{self.config.stage2.max_dnm_count}")
            samples_to_keep, n_removed = self._identify_outliers(data, warmup_params)
            logger.info(f"✓ Identified {n_removed} outlier individuals")
            self._log_memory_usage("After outlier identification")
        elif self.config.stage2.enabled and not warmup_params:
            logger.warning("Skipping outlier removal (no warmup parameters available)")
        else:
            logger.info("Stage 2 disabled in configuration")

        # Stage 3: Full optimisation with iterative outlier removal
        logger.info("")
        logger.info("="*60)
        logger.info("STAGE 3: FULL BAYESIAN OPTIMISATION")
        logger.info("="*60)
        final_params = {}
        scores = {}
        if self.config.stage3.enabled:
            logger.info(f"Running full optimisation with {self.config.stage3.n_trials} trials...")
            logger.info(f"Sampler: {self.config.stage3.sampler}, Pruner: {self.config.stage3.pruner}")
            
            # Iterative optimization: remove outliers based on current best params until stable
            max_iterations = 5
            for iteration in range(max_iterations):
                logger.info(f"Optimisation iteration {iteration + 1}...")
                
                final_params, scores = self._run_full_optimisation_sequential(
                    data, targets_by_type, samples_to_keep
                )
                
                # Check for new outliers under final params
                if self.config.stage2.enabled and final_params:
                    new_samples_to_keep, n_new_outliers = self._identify_outliers_with_params(
                        data, final_params, samples_to_keep
                    )
                    if n_new_outliers > 0:
                        logger.info(f"Found {n_new_outliers} new outliers under optimised parameters")
                        samples_to_keep = new_samples_to_keep
                        n_removed += n_new_outliers
                        logger.info(f"Removed {n_new_outliers} outliers, re-running optimisation...")
                    else:
                        logger.info(f"✓ No new outliers found - optimisation converged after {iteration + 1} iteration(s)")
                        break
                else:
                    break
            else:
                logger.warning(f"Optimisation did not converge after {max_iterations} iterations")
            
            logger.info(f"✓ Full optimisation complete")
            self._log_memory_usage("After full optimisation")
        else:
            logger.warning("Stage 3 disabled in configuration, using warmup parameters as final")
            final_params = warmup_params
            scores = {}

        # Create result object
        result = OptimisationResult(
            best_params=final_params,
            best_scores=scores,
            warmup_params=warmup_params if warmup_params else None,
            n_individuals_removed=n_removed,
            study=None
        )

        # Generate plots if requested
        if generate_plots and output_dir and final_params:
            logger.info("Generating optimization result plots")
            try:
                # Filter data for plotting
                if samples_to_keep is not None:
                    plot_df = data.variants[data.variants['SAMPLE'].isin(samples_to_keep)]
                else:
                    plot_df = data.variants
                
                plot_optimization_results(
                    data_df=plot_df,
                    reference_df=reference.variants,
                    best_params=final_params,
                    output_dir=output_dir,
                    save_filtered=True
                )
                save_parameters(final_params, output_dir)
            except Exception as e:
                logger.warning(f"Failed to generate plots: {e}")

        return result
    
    def _log_memory_usage(self, stage: str):
        """Log current memory usage."""
        try:
            import psutil
            process = psutil.Process()
            mem_mb = process.memory_info().rss / 1024 / 1024
            logger.info(f"Memory usage ({stage}): {mem_mb:.1f} MB")
        except ImportError:
            pass
    
    def _calculate_targets(self, reference: VariantDataset) -> Dict[str, np.ndarray]:
        """Calculate regression targets from reference data."""
        targets = {}
        
        for var_type in self.config.optimisation.variant_types:
            # Use boolean mask instead of copy
            type_mask = reference.variants['var_type'] == var_type
            if type_mask.sum() == 0:
                continue
            
            ref_subset = reference.variants[type_mask]
            
            # Count DNMs per sample
            counts = ref_subset.groupby('SAMPLE').size().rename('dnm_count').reset_index()

            # Get parental ages
            ages = ref_subset[['SAMPLE', 'paternal_age', 'maternal_age']].drop_duplicates()
            regression_data = ages.merge(counts, on='SAMPLE', how='left')
            regression_data = regression_data.dropna(subset=['paternal_age', 'maternal_age'])
            regression_data['dnm_count'] = regression_data['dnm_count'].fillna(0)
            
            try:
                model = smf.ols(self.config.optimisation.regression_formula, data=regression_data).fit()
                targets[var_type] = model.params.values
                logger.info(f"Calculated targets for {var_type}: {model.params.to_dict()}")
            except Exception as e:
                logger.warning(f"Failed to calculate targets for {var_type}: {e}")
        
        return targets
    
    def _get_filter_columns(self, var_type: str) -> List[str]:
        """Get list of columns that will be filtered for a variant type."""
        columns = []
        opt_columns = self.config.optimisation.get_optimisation_columns_for_variant_type(var_type)
        
        for col_config in opt_columns:
            columns.append(col_config.name)
        
        # Add linked columns
        linked_groups = self.config.optimisation.get_linked_column_groups()
        for group in linked_groups:
            for col_config in group:
                if col_config.name not in columns:
                    columns.append(col_config.name)
        
        return columns
    
    def _apply_range_prefilter(self, df: pd.DataFrame, var_type: str) -> pd.DataFrame:
        """Pre-filter DataFrame based on range constraints.
        
        This is particularly important for CMA-ES which treats all parameters as
        floats and cannot enforce hard range constraints during optimization.
        By pre-filtering, we ensure variants outside the valid range are excluded.
        """
        opt_columns = self.config.optimisation.get_optimisation_columns_for_variant_type(var_type)
        mask = np.ones(len(df), dtype=bool)
        
        for col_config in opt_columns:
            if col_config.range_constraint and col_config.optimisation == 'range':
                col = col_config.name
                
                # Handle computed columns (like abs_length)
                if col_config.computed and col in df.columns:
                    col_data = pd.to_numeric(df[col], errors='coerce')
                elif col in df.columns:
                    col_data = pd.to_numeric(df[col], errors='coerce')
                else:
                    continue
                
                # Apply range constraint
                range_min = col_config.range_constraint.min
                range_max = col_config.range_constraint.max
                
                col_mask = (col_data >= range_min) & (col_data <= range_max) & (~col_data.isna())
                n_filtered = (~col_mask).sum()
                
                if n_filtered > 0:
                    logger.info(f"  Pre-filtering {col}: removing {n_filtered} variants outside [{range_min}, {range_max}]")
                
                mask &= col_mask.values
        
        return df[mask]
    
    def _run_warmup_sequential(
        self,
        data: VariantDataset,
        targets: Dict[str, np.ndarray]
    ) -> Dict[str, Dict[str, Any]]:
        """Stage 1: Run warmup optimization for each variant type.
        
        Uses Optuna's n_jobs for parallel trial evaluation within each variant type.
        """
        warmup_params = {}
        is_cmaes = self.config.stage3.sampler == 'cmaes'
        
        for var_type in targets:
            logger.info(f"Processing {var_type}...")
            
            # Get filter columns for this variant type
            filter_columns = self._get_filter_columns(var_type)
            
            # Filter to this variant type using boolean mask
            type_mask = data.variants['var_type'] == var_type
            var_df = data.variants[type_mask].copy()
            
            if len(var_df) == 0:
                continue
            
            # For CMA-ES, pre-filter based on range constraints
            if is_cmaes:
                var_df = self._apply_range_prefilter(var_df, var_type)
                if len(var_df) == 0:
                    logger.warning(f"  No variants remaining after range pre-filter for {var_type}")
                    continue
            
            # Extract numpy arrays
            logger.info(f"  Extracting arrays for {len(var_df)} {var_type} variants...")
            arrays = extract_optimisation_arrays(var_df, filter_columns)
            
            # Run optimization with parallel trials
            params, _ = self._optimize_with_arrays(
                arrays,
                targets[var_type],
                self.config.stage1.n_trials,
                f"warmup_{var_type}",
                var_type
            )
            
            if params:
                warmup_params[var_type] = params
            
            # Explicit cleanup
            del arrays
            gc.collect()
            self._log_memory_usage(f"After {var_type} warmup")
        
        return warmup_params
    
    def _identify_outliers(
        self,
        data: VariantDataset,
        warmup_params: Dict[str, Dict[str, Any]]
    ) -> Tuple[set, int]:
        """Stage 2: Identify outlier samples based on filtered DNM counts.
        
        Returns set of samples to keep and count of removed samples.
        """
        logger.info("Counting filtered DNMs per individual...")
        
        # Count filtered DNMs per sample across all variant types
        sample_counts = {}
        
        for var_type, params in warmup_params.items():
            type_mask = data.variants['var_type'] == var_type
            var_df = data.variants[type_mask]
            
            if len(var_df) == 0:
                continue
            
            logger.info(f"  {var_type}: {len(var_df)} variants before warmup filters")
            
            # Apply filters using boolean mask
            filter_mask = np.ones(len(var_df), dtype=bool)
            for param, value in params.items():
                if param.startswith('min_'):
                    col = param[4:]
                    if col in var_df.columns:
                        col_vals = pd.to_numeric(var_df[col], errors='coerce')
                        filter_mask &= (col_vals >= value).values & (~col_vals.isna()).values
                elif param.startswith('max_'):
                    col = param[4:]
                    if col in var_df.columns:
                        col_vals = pd.to_numeric(var_df[col], errors='coerce')
                        filter_mask &= (col_vals <= value).values & (~col_vals.isna()).values
            
            filtered_count = filter_mask.sum()
            logger.info(f"  {var_type}: {filtered_count} variants after warmup filters ({filtered_count/len(var_df)*100:.1f}% retained)")
            
            # Count per sample
            filtered_df = var_df[filter_mask]
            type_counts = filtered_df.groupby('SAMPLE').size()
            
            for sample, count in type_counts.items():
                sample_counts[sample] = sample_counts.get(sample, 0) + count
        
        # Convert to series for analysis
        total_counts = pd.Series(sample_counts).sort_values()
        initial_count = len(total_counts)
        
        logger.info(f"DNM count distribution across {initial_count} individuals:")
        logger.info(f"  Min: {total_counts.min()}, Max: {total_counts.max()}, Mean: {total_counts.mean():.1f}, Median: {total_counts.median():.1f}")
        
        # Identify samples to keep
        samples_to_keep = set(total_counts.index)
        
        if self.config.stage2.min_dnm_count:
            below_min = total_counts[total_counts < self.config.stage2.min_dnm_count]
            samples_to_keep -= set(below_min.index)
            logger.info(f"Removing {len(below_min)} individuals with < {self.config.stage2.min_dnm_count} DNMs")
        
        if self.config.stage2.max_dnm_count:
            above_max = total_counts[total_counts > self.config.stage2.max_dnm_count]
            samples_to_keep -= set(above_max.index)
            logger.info(f"Removing {len(above_max)} individuals with > {self.config.stage2.max_dnm_count} DNMs")
        
        n_removed = initial_count - len(samples_to_keep)
        logger.info(f"Outlier removal complete: kept {len(samples_to_keep)}/{initial_count} individuals")
        
        return samples_to_keep, n_removed
    
    def _identify_outliers_with_params(
        self,
        data: VariantDataset,
        params: Dict[str, Dict[str, Any]],
        current_samples: Optional[set] = None
    ) -> Tuple[set, int]:
        """Identify outlier samples based on filtered DNM counts with given params.
        
        Similar to _identify_outliers but works with any parameter set and 
        optionally restricts to a subset of samples.
        
        Returns set of samples to keep and count of removed samples.
        """
        # Count filtered DNMs per sample across all variant types
        sample_counts = {}
        
        for var_type, var_params in params.items():
            type_mask = data.variants['var_type'] == var_type
            if current_samples is not None:
                type_mask &= data.variants['SAMPLE'].isin(current_samples)
            
            var_df = data.variants[type_mask]
            
            if len(var_df) == 0:
                continue
            
            # Apply filters using boolean mask
            filter_mask = np.ones(len(var_df), dtype=bool)
            for param, value in var_params.items():
                if param.startswith('min_'):
                    col = param[4:]
                    if col in var_df.columns:
                        col_vals = pd.to_numeric(var_df[col], errors='coerce')
                        filter_mask &= (col_vals >= value).values & (~col_vals.isna()).values
                elif param.startswith('max_'):
                    col = param[4:]
                    if col in var_df.columns:
                        col_vals = pd.to_numeric(var_df[col], errors='coerce')
                        filter_mask &= (col_vals <= value).values & (~col_vals.isna()).values
            
            # Count per sample
            filtered_df = var_df[filter_mask]
            type_counts = filtered_df.groupby('SAMPLE').size()
            
            for sample, count in type_counts.items():
                sample_counts[sample] = sample_counts.get(sample, 0) + count
        
        # Convert to series for analysis
        total_counts = pd.Series(sample_counts).sort_values()
        
        # Start with current samples or all samples
        if current_samples is not None:
            samples_to_keep = set(current_samples)
        else:
            samples_to_keep = set(total_counts.index)
        
        initial_count = len(samples_to_keep)
        
        # Apply thresholds
        if self.config.stage2.min_dnm_count:
            below_min = total_counts[total_counts < self.config.stage2.min_dnm_count]
            samples_to_keep -= set(below_min.index)
        
        if self.config.stage2.max_dnm_count:
            above_max = total_counts[total_counts > self.config.stage2.max_dnm_count]
            removed_above = set(above_max.index) & samples_to_keep
            if removed_above:
                logger.info(f"  Samples exceeding max_dnm_count ({self.config.stage2.max_dnm_count}): {len(removed_above)}")
                logger.info(f"  Their counts: {[int(total_counts[s]) for s in list(removed_above)[:10]]}{'...' if len(removed_above) > 10 else ''}")
            samples_to_keep -= set(above_max.index)
        
        n_removed = initial_count - len(samples_to_keep)
        
        return samples_to_keep, n_removed
    
    def _run_full_optimisation_sequential(
        self,
        data: VariantDataset,
        targets: Dict[str, np.ndarray],
        samples_to_keep: Optional[set] = None
    ) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, float]]:
        """Stage 3: Run full optimization for each variant type.
        
        Uses Optuna's n_jobs for parallel trial evaluation within each variant type.
        """
        final_params = {}
        scores = {}
        is_cmaes = self.config.stage3.sampler == 'cmaes'
        
        for var_type in targets:
            logger.info(f"Processing {var_type}...")
            
            # Get filter columns
            filter_columns = self._get_filter_columns(var_type)
            
            # Filter to this variant type and kept samples
            type_mask = data.variants['var_type'] == var_type
            if samples_to_keep is not None:
                type_mask &= data.variants['SAMPLE'].isin(samples_to_keep)
            
            var_df = data.variants[type_mask].copy()
            
            if len(var_df) == 0:
                continue
            
            # For CMA-ES, pre-filter based on range constraints
            if is_cmaes:
                var_df = self._apply_range_prefilter(var_df, var_type)
                if len(var_df) == 0:
                    logger.warning(f"  No variants remaining after range pre-filter for {var_type}")
                    continue
            
            # Extract numpy arrays
            logger.info(f"  Extracting arrays for {len(var_df)} {var_type} variants...")
            arrays = extract_optimisation_arrays(var_df, filter_columns)
            
            # Run optimization with parallel trials
            params, best_score = self._optimize_with_arrays(
                arrays,
                targets[var_type],
                self.config.stage3.n_trials,
                f"final_{var_type}",
                var_type,
                use_pruning=self.config.stage3.pruner is not None
            )
            
            if params:
                final_params[var_type] = params
                scores[var_type] = best_score
            
            # Explicit cleanup
            del arrays
            gc.collect()
            self._log_memory_usage(f"After {var_type} full optimisation")
        
        return final_params, scores
    
    def _optimize_with_arrays(
        self,
        arrays: OptimisationArrays,
        targets: np.ndarray,
        n_trials: int,
        study_name: str,
        var_type: str,
        use_pruning: bool = False
    ) -> Tuple[Optional[Dict[str, Any]], float]:
        """Optimize using pre-extracted numpy arrays (memory efficient).
        
        Uses Optuna's n_jobs for parallel trial evaluation via multiprocessing.
        """
        
        if arrays.n_variants == 0:
            return None, float('inf')
        
        # Pre-compute regression data structure (only samples with valid ages)
        valid_sample_indices = np.where(arrays.valid_age_mask)[0]
        n_valid_samples = len(valid_sample_indices)
        
        if n_valid_samples < 3:
            logger.warning(f"Not enough samples with valid ages: {n_valid_samples}")
            return None, float('inf')
        
        # Pre-extract ages for valid samples
        valid_paternal_ages = arrays.paternal_ages[valid_sample_indices]
        valid_maternal_ages = arrays.maternal_ages[valid_sample_indices]
        
        def objective(trial):
            # Suggest parameters
            params = self._suggest_params_for_arrays(trial, arrays, var_type)
            
            # Apply filters using boolean mask (no copy!)
            mask = apply_filters_mask(arrays, params)
            n_retained = mask.sum()
            
            if n_retained == 0:
                return float('inf')
            
            # Count per sample using fast numpy bincount
            counts = count_per_sample_fast(arrays.sample_ids, mask, arrays.n_samples)
            
            # Get counts for valid samples only
            valid_counts = counts[valid_sample_indices].astype(np.float64)
            
            # Check if we have enough variation
            if valid_counts.sum() == 0:
                return float('inf')
            
            try:
                # Build regression data (reuse pre-allocated arrays)
                regression_data = pd.DataFrame({
                    'dnm_count': valid_counts,
                    'paternal_age': valid_paternal_ages,
                    'maternal_age': valid_maternal_ages
                })
                
                # Fit regression
                model = smf.ols(self.config.optimisation.regression_formula, data=regression_data).fit()
                fitted_intercept = model.params.values[0]
                fitted_slope = model.params.values[1] if len(model.params.values) > 1 else 0
                target_intercept = targets[0]
                target_slope = targets[1] if len(targets) > 1 else 0
                
                # Calculate loss using relative error
                weights = self.config.optimisation.get_regression_weights(var_type)
                intercept_weight = weights[0] if len(weights) > 0 else 1.0
                slope_weight = weights[1] if len(weights) > 1 else 1.0
                
                slope_rel_error = (fitted_slope - target_slope) / (np.abs(target_slope) + 1e-10)
                intercept_rel_error = (fitted_intercept - target_intercept) / (np.abs(target_intercept) + 1e-10)
                
                total_loss = (intercept_weight * intercept_rel_error ** 2) + (slope_weight * slope_rel_error ** 2)
                
                if use_pruning and trial.number > 0:
                    trial.report(total_loss, trial.number)
                    if trial.should_prune():
                        raise optuna.TrialPruned()
                
                return total_loss
                
            except Exception as e:
                logger.debug(f"Trial {trial.number}: Regression failed: {e}")
                return float('inf')
        
        # Create study
        sampler_map = {
            'tpe': optuna.samplers.TPESampler(seed=self.config.seed, multivariate=self.config.stage3.multivariate),
            'cmaes': optuna.samplers.CmaEsSampler(seed=self.config.seed),
            'random': optuna.samplers.RandomSampler(seed=self.config.seed)
        }
        
        pruner_map = {
            'successive_halving': optuna.pruners.SuccessiveHalvingPruner(min_resource=1, reduction_factor=3),
            'hyperband': optuna.pruners.HyperbandPruner(),
            'median': optuna.pruners.MedianPruner()
        }
        
        sampler = sampler_map.get(self.config.stage3.sampler, optuna.samplers.TPESampler(seed=self.config.seed))
        pruner = pruner_map.get(self.config.stage3.pruner) if use_pruning else None
        
        study = optuna.create_study(
            direction='minimize',
            sampler=sampler,
            pruner=pruner,
            study_name=study_name
        )
        
        # Run optimisation with parallel trials
        # Note: n_jobs > 1 uses joblib multiprocessing which requires picklable objectives
        # The loky backend (joblib default) handles closures via cloudpickle
        n_jobs = self.config.max_workers
        if n_jobs > 1:
            logger.info(f"  Running {n_trials} trials with {n_jobs} parallel workers...")
        study.optimize(
            objective, 
            n_trials=n_trials, 
            n_jobs=n_jobs, 
            show_progress_bar=(n_jobs == 1)  # Progress bar doesn't work well with multiprocessing
        )
        
        score_str = f"{study.best_value:.4e}" if study.best_value < 0.01 else f"{study.best_value:.6f}"
        logger.info(f"Optimisation for {study_name} finished. Best score: {score_str}")
        
        # Log best params details
        self._log_best_params(study, arrays)
        
        return study.best_params, study.best_value
    
    def _log_best_params(self, study: optuna.Study, arrays: OptimisationArrays):
        """Log details about the best parameters found."""
        try:
            best_params = study.best_params
            mask = apply_filters_mask(arrays, best_params)
            n_retained = mask.sum()
            retention_rate = n_retained / arrays.n_variants * 100
            
            logger.info(f"  Retention: {n_retained}/{arrays.n_variants} variants ({retention_rate:.1f}%)")
            logger.info(f"  Best params: {best_params}")
        except Exception as e:
            logger.debug(f"Could not log best params: {e}")
    
    def _suggest_params_for_arrays(
        self,
        trial: optuna.Trial,
        arrays: OptimisationArrays,
        var_type: str
    ) -> Dict[str, Any]:
        """Suggest parameters using pre-extracted arrays."""
        params = {}
        
        opt_columns = self.config.optimisation.get_optimisation_columns_for_variant_type(var_type)
        if not opt_columns:
            return params
        
        is_cmaes = self.config.stage3.sampler == 'cmaes'
        
        # Build skip columns for linked groups
        linked_groups = self.config.optimisation.get_linked_column_groups()
        skip_columns = set()
        for group in linked_groups:
            if len(group) >= 2:
                for col_config in group[1:]:
                    skip_columns.add(col_config.name)
        
        for col_config in opt_columns:
            col = col_config.name
            
            if col in skip_columns:
                continue
            
            # Check if column exists in arrays
            if not col_config.computed and col not in arrays.filter_arrays:
                continue
            
            # Get bounds
            if col_config.range_constraint:
                lower = col_config.range_constraint.min
                upper = col_config.range_constraint.max
            elif col_config.computed:
                continue
            else:
                col_data = arrays.filter_arrays[col]
                col_data = col_data[~np.isnan(col_data)]
                if len(col_data) == 0:
                    continue
                
                if col_config.optimisation == 'minimum':
                    lower = float(np.nanpercentile(col_data, 10))
                    upper = float(np.nanmax(col_data))
                else:  # maximum
                    lower = float(np.nanmin(col_data))
                    upper = float(np.nanpercentile(col_data, 90))
            
            if trial.number == 0:
                prefix = 'max' if col_config.optimisation == 'minimum' else 'min'
                logger.info(f"  {prefix}_{col} search space: [{lower:.2f}, {upper:.2f}]")
            
            # Suggest value
            if col_config.optimisation == 'minimum':
                if col_config.dtype == 'int' and not is_cmaes:
                    params[f'max_{col}'] = trial.suggest_int(f'max_{col}', int(lower), int(upper))
                else:
                    value = trial.suggest_float(f'max_{col}', float(lower), float(upper))
                    # Clip to bounds for CMA-ES (it can sample outside bounds)
                    value = max(lower, min(upper, value))
                    params[f'max_{col}'] = int(round(value)) if col_config.dtype == 'int' else value
            
            elif col_config.optimisation == 'maximum':
                if col_config.dtype == 'int' and not is_cmaes:
                    params[f'min_{col}'] = trial.suggest_int(f'min_{col}', int(lower), int(upper))
                else:
                    value = trial.suggest_float(f'min_{col}', float(lower), float(upper))
                    # Clip to bounds for CMA-ES
                    value = max(lower, min(upper, value))
                    params[f'min_{col}'] = int(round(value)) if col_config.dtype == 'int' else value
            
            elif col_config.optimisation == 'range':
                if hasattr(col_config.range_constraint, 'lower') and hasattr(col_config.range_constraint, 'scale'):
                    # Symmetric range
                    lower_bound = col_config.range_constraint.lower
                    upper_bound = col_config.range_constraint.scale - col_config.range_constraint.lower
                    
                    if col_config.dtype == 'int' and not is_cmaes:
                        suggested_lower = trial.suggest_int(f'min_{col}', int(lower_bound), int(upper_bound))
                    else:
                        value = trial.suggest_float(f'min_{col}', float(lower_bound), float(upper_bound))
                        # Clip to bounds for CMA-ES
                        value = max(lower_bound, min(upper_bound, value))
                        suggested_lower = int(round(value)) if col_config.dtype == 'int' else value
                    
                    suggested_upper = col_config.range_constraint.scale - suggested_lower
                    params[f'min_{col}'] = suggested_lower
                    params[f'max_{col}'] = suggested_upper
                else:
                    # Regular range
                    if col_config.dtype == 'int' and not is_cmaes:
                        params[f'min_{col}'] = trial.suggest_int(f'min_{col}', int(lower), int(upper))
                        params[f'max_{col}'] = trial.suggest_int(f'max_{col}', int(lower), int(upper))
                    else:
                        value_min = trial.suggest_float(f'min_{col}', float(lower), float(upper))
                        value_max = trial.suggest_float(f'max_{col}', float(lower), float(upper))
                        # Clip to bounds for CMA-ES
                        value_min = max(lower, min(upper, value_min))
                        value_max = max(lower, min(upper, value_max))
                        if col_config.dtype == 'int':
                            params[f'min_{col}'] = int(round(value_min))
                            params[f'max_{col}'] = int(round(value_max))
                        else:
                            params[f'min_{col}'] = value_min
                            params[f'max_{col}'] = value_max
        
        # Handle linked columns
        for group in linked_groups:
            if len(group) < 2:
                continue
            
            first_col = group[0].name
            param_keys = [k for k in params.keys() if k.endswith(f'_{first_col}')]
            
            for key in param_keys:
                param_value = params[key]
                prefix = key.split('_')[0]
                
                for linked_col in group[1:]:
                    linked_key = f'{prefix}_{linked_col.name}'
                    params[linked_key] = param_value
        
        return params


# Alias for backward compatibility
OptimisationPipeline = MemoryEfficientPipeline
