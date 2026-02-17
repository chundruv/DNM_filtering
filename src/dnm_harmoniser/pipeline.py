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
from .plotting import (
    plot_optimization_results, 
    save_parameters,
    plot_optuna_diagnostics,
    plot_regression_trajectory,
    plot_optimization_summary,
)
from .vaf_quality import compute_vaf_penalty


logger = logging.getLogger(__name__)


@dataclass
class OptimisationResult:
    """Container for optimisation results."""
    best_params: Dict[str, Dict[str, Any]]
    best_scores: Dict[str, float]
    warmup_params: Optional[Dict[str, Dict[str, Any]]]
    n_individuals_removed: int
    studies: Optional[Dict[str, optuna.Study]]  # Studies per variant type
    trial_history: Optional[Dict[str, List[Dict[str, Any]]]]  # Trial history per variant type
    targets: Optional[Dict[str, Dict[str, float]]]  # Target slope/intercept per variant type
    
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
    
    # VAF values per variant (for VAF quality penalty in objective)
    vaf_values: Optional[np.ndarray]  # shape (n_variants,), None if not available
    
    # Parental ages per sample (indexed by sample_idx)
    paternal_ages: np.ndarray  # shape (n_samples,)
    valid_age_mask: np.ndarray  # bool array for samples with valid ages
    
    # Original number of variants
    n_variants: int


def extract_optimisation_arrays(
    df: pd.DataFrame,
    filter_columns: List[str],
    sample_col: str = 'SAMPLE',
    all_sample_ages: Optional[pd.DataFrame] = None,
    vaf_col: Optional[str] = None,
) -> OptimisationArrays:
    """Extract numpy arrays from DataFrame for memory-efficient optimization.
    
    This function extracts only the data needed for the objective function
    as compact numpy arrays, allowing the original DataFrame to be freed.
    
    Args:
        df: DataFrame containing variants (may be empty if no variants of this type)
        filter_columns: List of column names to extract for filtering
        sample_col: Name of sample ID column
        all_sample_ages: Optional DataFrame with columns [SAMPLE, paternal_age]
                        containing ALL samples to include in regression (even those with
                        0 variants of this type). If None, samples are derived from df.
        vaf_col: Optional column name for VAF values (for VAF quality penalty)
    """
    # Create sample ID mapping - use all_sample_ages if provided, otherwise derive from df
    if all_sample_ages is not None:
        unique_samples = all_sample_ages['SAMPLE'].unique()
    else:
        unique_samples = df[sample_col].unique()
    
    sample_id_to_idx = {s: i for i, s in enumerate(unique_samples)}
    idx_to_sample_id = {i: s for s, i in sample_id_to_idx.items()}
    n_samples = len(unique_samples)
    
    # Encode sample IDs as integers (only for variants in df)
    # Samples not in df will have 0 count via bincount
    if len(df) > 0:
        # Map sample IDs, using -1 for any samples in df but not in our sample list
        # (shouldn't happen if all_sample_ages is comprehensive, but be safe)
        sample_ids = df[sample_col].map(lambda x: sample_id_to_idx.get(x, -1)).values.astype(np.int32)
        # Filter out any -1 values (samples not in our list)
        valid_variant_mask = sample_ids >= 0
        sample_ids = sample_ids[valid_variant_mask]
    else:
        sample_ids = np.array([], dtype=np.int32)
        valid_variant_mask = np.array([], dtype=bool)
    
    # Extract filter columns (only for valid variants)
    filter_arrays = {}
    for col in filter_columns:
        if col in df.columns:
            col_data = pd.to_numeric(df[col], errors='coerce').values.astype(np.float32)
            if len(valid_variant_mask) > 0:
                filter_arrays[col] = col_data[valid_variant_mask]
            else:
                filter_arrays[col] = col_data
    
    # Extract parental ages per sample
    if all_sample_ages is not None:
        # Use provided ages
        sample_ages = all_sample_ages.set_index('SAMPLE')
    else:
        # Derive from df
        sample_ages = df[[sample_col, 'paternal_age']].drop_duplicates(subset=[sample_col])
        sample_ages = sample_ages.set_index(sample_col)
    
    paternal_ages = np.full(n_samples, np.nan, dtype=np.float32)
    
    for sample, idx in sample_id_to_idx.items():
        if sample in sample_ages.index:
            paternal_ages[idx] = sample_ages.loc[sample, 'paternal_age']
    
    # Only require valid paternal age (maternal age not used in optimization)
    valid_age_mask = ~np.isnan(paternal_ages)
    
    # n_variants is the number of variants we're actually using (after filtering to known samples)
    n_variants = len(sample_ids)
    
    # Extract VAF values if column is specified and exists
    vaf_values = None
    if vaf_col and vaf_col in df.columns and len(df) > 0:
        vaf_data = pd.to_numeric(df[vaf_col], errors='coerce').values.astype(np.float32)
        if len(valid_variant_mask) > 0:
            vaf_values = vaf_data[valid_variant_mask]
        else:
            vaf_values = vaf_data
        logger.debug(f"Extracted {len(vaf_values)} VAF values from column '{vaf_col}'")
    
    return OptimisationArrays(
        sample_ids=sample_ids,
        sample_id_to_idx=sample_id_to_idx,
        idx_to_sample_id=idx_to_sample_id,
        n_samples=n_samples,
        filter_arrays=filter_arrays,
        vaf_values=vaf_values,
        paternal_ages=paternal_ages,
        valid_age_mask=valid_age_mask,
        n_variants=n_variants
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
            studies=self._studies if hasattr(self, '_studies') else None,
            trial_history=self._trial_history if hasattr(self, '_trial_history') else None,
            targets=self._targets_dict if hasattr(self, '_targets_dict') else None,
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
                    save_filtered=True,
                    config=self.config
                )
                save_parameters(final_params, output_dir)
                
                # Generate Optuna diagnostic plots
                if hasattr(self, '_studies') and self._studies:
                    try:
                        plot_optuna_diagnostics(self._studies, output_dir)
                    except Exception as e:
                        logger.warning(f"Failed to generate Optuna diagnostics: {e}")
                
                # Generate regression trajectory plots
                if hasattr(self, '_trial_history') and self._trial_history and hasattr(self, '_targets_dict'):
                    try:
                        plot_regression_trajectory(self._trial_history, self._targets_dict, output_dir)
                        plot_optimization_summary(self._studies, self._trial_history, self._targets_dict, output_dir)
                    except Exception as e:
                        logger.warning(f"Failed to generate trajectory plots: {e}")
                        
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
        """Calculate regression targets from reference data.
        
        Uses paternal_age only regression: dnm_count ~ paternal_age
        Returns coefficients as [Intercept, paternal_age_slope]
        
        All samples with valid paternal ages are included in the regression,
        even if they have 0 variants of a specific type.
        """
        targets = {}
        
        # Get ALL samples with valid paternal ages from reference (across all variant types)
        all_ref_ages = reference.variants[['SAMPLE', 'paternal_age']].drop_duplicates(subset=['SAMPLE'])
        all_ref_ages = all_ref_ages.dropna(subset=['paternal_age'])
        logger.info(f"Reference dataset has {len(all_ref_ages)} samples with valid paternal ages")
        
        for var_type in self.config.optimisation.variant_types:
            # Use boolean mask instead of copy
            type_mask = reference.variants['var_type'] == var_type
            ref_subset = reference.variants[type_mask]
            
            # Count DNMs per sample (only samples with variants will appear)
            if len(ref_subset) > 0:
                counts = ref_subset.groupby('SAMPLE').size().rename('dnm_count').reset_index()
            else:
                # No variants of this type - all samples have 0 count
                counts = pd.DataFrame({'SAMPLE': [], 'dnm_count': []})

            # Merge with ALL samples (to include those with 0 variants of this type)
            regression_data = all_ref_ages.merge(counts, on='SAMPLE', how='left')
            regression_data['dnm_count'] = regression_data['dnm_count'].fillna(0)
            
            n_with_variants = (regression_data['dnm_count'] > 0).sum()
            logger.info(f"Reference {var_type}: {n_with_variants}/{len(regression_data)} samples have >=1 variant")
            
            try:
                # Use paternal_age only formula
                model = smf.ols('dnm_count ~ paternal_age', data=regression_data).fit()
                targets[var_type] = model.params.values  # [Intercept, paternal_age]
                logger.info(f"Reference targets for {var_type}: intercept={model.params['Intercept']:.4f}, slope={model.params['paternal_age']:.4f}")
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
        All samples with valid ages are included in each variant type's regression,
        even if they have 0 variants of that type.
        """
        warmup_params = {}
        is_cmaes = self.config.stage3.sampler == 'cmaes'
        
        # Get ALL sample ages from the full dataset (no outlier removal yet in warmup)
        all_sample_ages = data.variants[['SAMPLE', 'paternal_age']].drop_duplicates(subset=['SAMPLE'])
        all_sample_ages = all_sample_ages.dropna(subset=['paternal_age'])  # Must have valid paternal age
        logger.info(f"Warmup: Total samples for regression (across all types): {len(all_sample_ages)}")
        
        for var_type in targets:
            logger.info(f"Processing {var_type}...")
            
            # Get filter columns for this variant type
            filter_columns = self._get_filter_columns(var_type)
            
            # Filter to this variant type using boolean mask
            type_mask = data.variants['var_type'] == var_type
            var_df = data.variants[type_mask].copy()
            
            # Note: var_df may be empty if no variants of this type exist,
            # but we still want to include all samples in regression (with 0 counts)
            
            # For CMA-ES, pre-filter based on range constraints
            if is_cmaes and len(var_df) > 0:
                var_df = self._apply_range_prefilter(var_df, var_type)
            
            # Extract numpy arrays - pass all_sample_ages to include samples with 0 variants
            n_variants = len(var_df)
            logger.info(f"  Extracting arrays for {n_variants} {var_type} variants ({len(all_sample_ages)} samples)...")
            # Pass vaf_col if VAF quality metric is configured
            vaf_col = self.config.optimisation.vaf_column if self.config.optimisation.vaf_quality_metric else None
            arrays = extract_optimisation_arrays(var_df, filter_columns, all_sample_ages=all_sample_ages, vaf_col=vaf_col)
            
            # Debug: Log how many samples have variants vs 0 BEFORE any filtering
            if arrays.n_variants > 0:
                all_counts = count_per_sample_fast(arrays.sample_ids, np.ones(arrays.n_variants, dtype=bool), arrays.n_samples)
            else:
                all_counts = np.zeros(arrays.n_samples, dtype=int)
            valid_idx = np.where(arrays.valid_age_mask)[0]
            valid_counts_pre = all_counts[valid_idx]
            n_with_var = np.sum(valid_counts_pre > 0)
            n_zero = np.sum(valid_counts_pre == 0)
            mean_count = np.mean(valid_counts_pre)
            logger.info(f"  DEBUG {var_type} Stage3: {n_with_var}/{len(valid_counts_pre)} have >=1 variant, {n_zero} have 0, mean={mean_count:.2f}")
            
            # Debug: Log how many samples have variants vs 0
            all_counts = count_per_sample_fast(arrays.sample_ids, np.ones(arrays.n_variants, dtype=bool), arrays.n_samples)
            valid_idx = np.where(arrays.valid_age_mask)[0]
            valid_counts = all_counts[valid_idx]
            n_with_var = np.sum(valid_counts > 0)
            n_zero = np.sum(valid_counts == 0)
            logger.info(f"  DEBUG {var_type}: Before filtering - {n_with_var}/{len(valid_counts)} have >=1 variant, {n_zero} have 0")
            
            # Run optimization with parallel trials
            params, _, _, _ = self._optimize_with_arrays(
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
        All non-outlier samples are included in each variant type's regression,
        even if they have 0 variants of that type.
        """
        final_params = {}
        scores = {}
        is_cmaes = self.config.stage3.sampler == 'cmaes'
        
        # Initialize storage for studies and trial history
        self._studies = {}
        self._trial_history = {}
        self._targets_dict = {}
        
        # Convert targets to dict format for plotting
        for var_type, target_arr in targets.items():
            self._targets_dict[var_type] = {
                'intercept': target_arr[0],
                'slope': target_arr[1] if len(target_arr) > 1 else 0
            }
        
        # Get ALL sample ages from the full dataset (filtered to non-outliers)
        # This ensures samples with 0 variants for a specific type are still included
        if samples_to_keep is not None:
            all_samples_df = data.variants[data.variants['SAMPLE'].isin(samples_to_keep)]
        else:
            all_samples_df = data.variants
        
        all_sample_ages = all_samples_df[['SAMPLE', 'paternal_age']].drop_duplicates(subset=['SAMPLE'])
        all_sample_ages = all_sample_ages.dropna(subset=['paternal_age'])  # Must have valid paternal age
        logger.info(f"Total samples for regression (across all types): {len(all_sample_ages)}")
        
        for var_type in targets:
            logger.info(f"Processing {var_type}...")
            
            # Get filter columns
            filter_columns = self._get_filter_columns(var_type)
            
            # Filter to this variant type and kept samples
            type_mask = data.variants['var_type'] == var_type
            if samples_to_keep is not None:
                type_mask &= data.variants['SAMPLE'].isin(samples_to_keep)
            
            var_df = data.variants[type_mask].copy()
            
            # Note: var_df may be empty if no variants of this type exist,
            # but we still want to include all samples in regression (with 0 counts)
            
            # For CMA-ES, pre-filter based on range constraints
            if is_cmaes and len(var_df) > 0:
                var_df = self._apply_range_prefilter(var_df, var_type)
            
            # Extract numpy arrays - pass all_sample_ages to include samples with 0 variants
            n_variants = len(var_df)
            logger.info(f"  Extracting arrays for {n_variants} {var_type} variants ({len(all_sample_ages)} samples)...")
            # Pass vaf_col if VAF quality metric is configured
            vaf_col = self.config.optimisation.vaf_column if self.config.optimisation.vaf_quality_metric else None
            arrays = extract_optimisation_arrays(var_df, filter_columns, all_sample_ages=all_sample_ages, vaf_col=vaf_col)
            
            # Run optimization with parallel trials
            params, best_score, study, trial_hist = self._optimize_with_arrays(
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
                self._studies[var_type] = study
                self._trial_history[var_type] = trial_hist
            
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
    ) -> Tuple[Optional[Dict[str, Any]], float, Optional[optuna.Study], List[Dict[str, Any]]]:
        """
        Optimize filter parameters to match reference regression coefficients.
        
        This method directly optimizes the slope and intercept of the 
        dnm_count ~ paternal_age regression to match reference values.
        
        Returns:
            Tuple of (best_params, best_score, study, trial_history)
        """
        
        trial_history = []
        
        if arrays.n_variants == 0:
            return None, float('inf'), None, trial_history
        
        # 1. Setup: Identify samples with valid ages
        valid_sample_indices = np.where(arrays.valid_age_mask)[0]
        n_valid_samples = len(valid_sample_indices)
        
        if n_valid_samples < 10:
            logger.warning(f"Not enough samples with valid ages: {n_valid_samples}")
            return None, float('inf'), None, trial_history
        
        # 2. Extract Reference Coefficients
        # Assumes targets are [Intercept, Paternal_Slope, ...]
        ref_intercept = targets[0]
        ref_pat_slope = targets[1] if len(targets) > 1 else 0
        
        # Get weights from config [intercept_weight, slope_weight]
        weights = self.config.optimisation.get_regression_weights(var_type)
        intercept_weight = weights[0] if len(weights) > 0 else 1.0
        slope_weight = weights[1] if len(weights) > 1 else 1.0
        
        # Get loss function from config
        loss_fn = self.config.optimisation.get_loss_function(var_type)
        huber_delta = self.config.optimisation.huber_delta
        asymmetric_penalty = self.config.optimisation.asymmetric_penalty
        intercept_tolerance = self.config.optimisation.intercept_tolerance
        
        logger.info(f"  Reference targets: intercept={ref_intercept:.4f}, slope={ref_pat_slope:.4f}")
        logger.info(f"  === LOSS CONFIG ===")
        logger.info(f"      Loss function: {loss_fn}")
        logger.info(f"      Weights: intercept={intercept_weight:.2f}, slope={slope_weight:.2f}")
        if loss_fn == "huber":
            logger.info(f"      Huber delta: {huber_delta}")
        elif loss_fn == "asymmetric":
            logger.info(f"      Asymmetric penalty: {asymmetric_penalty}x")
        elif loss_fn in ("slope_priority", "max_slope"):
            logger.info(f"      Intercept tolerance: {intercept_tolerance:.1%}")
        logger.info(f"  ====================")
        
        # VAF quality penalty configuration
        vaf_metric = self.config.optimisation.vaf_quality_metric
        vaf_weight = self.config.optimisation.get_vaf_quality_weight(var_type)
        vaf_min = self.config.optimisation.vaf_min
        vaf_max = self.config.optimisation.vaf_max
        has_vaf_penalty = (vaf_metric is not None and vaf_weight > 0 
                          and arrays.vaf_values is not None)
        if has_vaf_penalty:
            logger.info(f"  VAF quality penalty: metric={vaf_metric}, weight={vaf_weight:.3f}")
        else:
            if vaf_metric is not None and arrays.vaf_values is None:
                logger.warning(f"  VAF quality metric '{vaf_metric}' configured but no VAF values available")

        # 3. Pre-extract paternal ages for valid samples
        valid_pat_ages = arrays.paternal_ages[valid_sample_indices]
        
        # Pre-compute values for fast OLS
        pat_mean = np.mean(valid_pat_ages)
        pat_centered = valid_pat_ages - pat_mean
        pat_var_sum = np.sum(pat_centered ** 2)

        def objective(trial):
            # A. Suggest Parameters
            params = self._suggest_params_for_arrays(trial, arrays, var_type)
            
            # B. Apply Filters
            mask = apply_filters_mask(arrays, params)
            
            # C. Count variants per sample
            counts = count_per_sample_fast(arrays.sample_ids, mask, arrays.n_samples)
            valid_counts = counts[valid_sample_indices].astype(np.float64)
            
            # D. Check for degenerate solutions
            if valid_counts.sum() == 0:
                return 1e9
            
            # E. Fast OLS: dnm_count ~ paternal_age
            count_mean = np.mean(valid_counts)
            count_centered = valid_counts - count_mean
            
            # slope = cov(x,y) / var(x)
            slope = np.sum(pat_centered * count_centered) / pat_var_sum
            intercept = count_mean - slope * pat_mean
            
            # F. Calculate loss based on selected loss function
            # Normalized errors (used by most loss functions)
            slope_rel_error = (slope - ref_pat_slope) / ref_pat_slope if ref_pat_slope != 0 else slope
            intercept_rel_error = (intercept - ref_intercept) / ref_intercept if ref_intercept != 0 else intercept
            
            if loss_fn == "weighted_mse":
                # Default: weighted mean squared error of normalized values
                slope_error = slope_rel_error ** 2
                intercept_error = intercept_rel_error ** 2
                loss = intercept_weight * intercept_error + slope_weight * slope_error
                
            elif loss_fn == "absolute":
                # L1 loss - less sensitive to large deviations
                slope_error = abs(slope_rel_error)
                intercept_error = abs(intercept_rel_error)
                loss = intercept_weight * intercept_error + slope_weight * slope_error
                
            elif loss_fn == "huber":
                # Huber loss - quadratic for small errors, linear for large
                def huber(x, delta):
                    if abs(x) <= delta:
                        return 0.5 * x ** 2
                    else:
                        return delta * (abs(x) - 0.5 * delta)
                slope_error = huber(slope_rel_error, huber_delta)
                intercept_error = huber(intercept_rel_error, huber_delta)
                loss = intercept_weight * intercept_error + slope_weight * slope_error
                
            elif loss_fn == "log_ratio":
                # Log ratio - good when values span orders of magnitude
                # Avoid log of negative/zero
                if slope > 0 and ref_pat_slope > 0:
                    slope_error = (np.log(slope) - np.log(ref_pat_slope)) ** 2
                else:
                    slope_error = slope_rel_error ** 2
                if intercept > 0 and ref_intercept > 0:
                    intercept_error = (np.log(intercept) - np.log(ref_intercept)) ** 2
                else:
                    intercept_error = intercept_rel_error ** 2
                loss = intercept_weight * intercept_error + slope_weight * slope_error
                
            elif loss_fn == "slope_priority":
                # Optimize slope first, only penalize intercept if way off
                slope_error = slope_rel_error ** 2
                # Soft constraint on intercept - only penalize if outside tolerance
                intercept_deviation = abs(intercept_rel_error)
                if intercept_deviation > intercept_tolerance:
                    intercept_error = (intercept_deviation - intercept_tolerance) ** 2
                else:
                    intercept_error = 0
                loss = slope_weight * slope_error + intercept_weight * intercept_error
                
            elif loss_fn == "asymmetric":
                # Penalize undershooting slope more than overshooting
                slope_error = slope_rel_error ** 2
                if slope < ref_pat_slope:
                    # Undershooting - apply extra penalty
                    slope_error *= asymmetric_penalty
                intercept_error = intercept_rel_error ** 2
                loss = intercept_weight * intercept_error + slope_weight * slope_error
                
            elif loss_fn == "max_slope":
                # Maximize slope with soft intercept constraint
                # Negative slope means we want to maximize it
                slope_loss = -slope / ref_pat_slope  # Negative to maximize
                
                # Add penalty if intercept deviates too much
                intercept_deviation = abs(intercept_rel_error)
                if intercept_deviation > intercept_tolerance:
                    intercept_penalty = intercept_weight * (intercept_deviation - intercept_tolerance) ** 2
                else:
                    intercept_penalty = 0
                
                loss = slope_weight * slope_loss + intercept_penalty
                # For tracking, use absolute values
                slope_error = slope_rel_error ** 2
                intercept_error = intercept_rel_error ** 2
                
            elif loss_fn == "correlation":
                # Maximize correlation between filtered counts and paternal age
                # This directly measures how well the paternal age effect is preserved
                if np.std(valid_counts) > 0:
                    corr = np.corrcoef(valid_counts, valid_pat_ages)[0, 1]
                    if np.isnan(corr):
                        corr = 0
                else:
                    corr = 0
                
                # Transform to loss: (1 - corr) so perfect correlation = 0, no correlation = 1
                # This keeps loss positive and intuitive
                loss = 1 - corr
                
                # For tracking
                slope_error = slope_rel_error ** 2
                intercept_error = intercept_rel_error ** 2
                
            else:
                # Fallback to weighted_mse
                slope_error = slope_rel_error ** 2
                intercept_error = intercept_rel_error ** 2
                loss = intercept_weight * intercept_error + slope_weight * slope_error
            
            # G. VAF quality penalty (if configured)
            vaf_penalty = 0.0
            if has_vaf_penalty:
                filtered_vaf = arrays.vaf_values[mask]
                if len(filtered_vaf) >= 10:
                    vaf_penalty = compute_vaf_penalty(
                        filtered_vaf, var_type, metric=vaf_metric,
                        vaf_min=vaf_min, vaf_max=vaf_max
                    )
                    loss = loss + vaf_weight * vaf_penalty
            
            # H. Track trial history for plotting
            trial_history.append({
                'trial_number': trial.number,
                'slope': slope,
                'intercept': intercept,
                'loss': loss,
                'slope_error': slope_error,
                'intercept_error': intercept_error,
                'vaf_penalty': vaf_penalty,
            })
            
            return loss
        
        # Create study
        sampler_map = {
            'tpe': optuna.samplers.TPESampler(seed=self.config.seed, multivariate=True), # Multivariate helps TPE find correlations
            'cmaes': optuna.samplers.CmaEsSampler(seed=self.config.seed),
            'random': optuna.samplers.RandomSampler(seed=self.config.seed)
        }
        
        sampler = sampler_map.get(self.config.stage3.sampler, optuna.samplers.TPESampler(seed=self.config.seed))
        
        # Disable pruner for stability unless explicitly requested and safe
        pruner = None 
        
        study = optuna.create_study(
            direction='minimize',
            sampler=sampler,
            pruner=pruner,
            study_name=study_name
        )
        
        n_jobs = self.config.max_workers
        # Cap workers to avoid overhead on simple objective
        if n_jobs > 1:
            logger.info(f"  Running {n_trials} trials...")
            
        study.optimize(
            objective, 
            n_trials=n_trials, 
            n_jobs=n_jobs, 
            show_progress_bar=(n_jobs == 1)
        )
        
        # Logging
        best_loss = study.best_value
        logger.info(f"Optimisation finished. Best loss: {best_loss:.6f}")
        
        # Log achieved slope/intercept
        best_params = study.best_params
        mask = apply_filters_mask(arrays, best_params)
        counts = count_per_sample_fast(arrays.sample_ids, mask, arrays.n_samples)
        valid_counts = counts[valid_sample_indices].astype(np.float64)
        
        # Detailed logging for debugging mismatch with plots
        n_with_variants = np.sum(valid_counts > 0)
        n_zero = np.sum(valid_counts == 0)
        n_total = len(valid_counts)
        mean_count = np.mean(valid_counts)
        logger.info(f"  OPTIMIZATION DATA: {n_with_variants}/{n_total} samples have >=1 variant, {n_zero} have 0, mean={mean_count:.2f}")
        
        count_mean = np.mean(valid_counts)
        count_centered = valid_counts - count_mean
        achieved_slope = np.sum(pat_centered * count_centered) / pat_var_sum
        achieved_intercept = count_mean - achieved_slope * pat_mean
        
        logger.info(f"  OPTIMIZATION REGRESSION: intercept={achieved_intercept:.4f} (target={ref_intercept:.4f}), slope={achieved_slope:.4f} (target={ref_pat_slope:.4f})")
        
        # Log VAF quality penalty for best params
        if has_vaf_penalty:
            filtered_vaf = arrays.vaf_values[mask]
            if len(filtered_vaf) >= 10:
                best_vaf_penalty = compute_vaf_penalty(
                    filtered_vaf, var_type, metric=vaf_metric,
                    vaf_min=vaf_min, vaf_max=vaf_max
                )
                logger.info(f"  VAF QUALITY: metric={vaf_metric}, penalty={best_vaf_penalty:.4f}, "
                           f"weighted_contribution={vaf_weight * best_vaf_penalty:.4f}")
        
        self._log_best_params(study, arrays)
        
        return study.best_params, study.best_value, study, trial_history
    
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