"""Optimized three-stage pipeline implementation."""

import pandas as pd
import numpy as np
import optuna
import statsmodels.formula.api as smf
from typing import Dict, List, Any, Optional, Tuple
from joblib import Memory
from concurrent.futures import ProcessPoolExecutor
import hashlib
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
            lines.append(f"  Best score: {self.best_scores.get(var_type, 0):.4f}")
            for param, value in params.items():
                if isinstance(value, float):
                    lines.append(f"  {param}: {value:.4f}")
                else:
                    lines.append(f"  {param}: {value}")
        if self.n_individuals_removed > 0:
            lines.append(f"\nOutliers removed: {self.n_individuals_removed}")
        return '\n'.join(lines)


class OptimisationPipeline:
    """Three-stage optimisation pipeline with caching and parallelization."""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        
        # Set up caching
        if config.use_cache:
            cache_dir = config.cache_dir or './cache'
            self.memory = Memory(cache_dir, verbose=0)
            self._stage1_cached = self.memory.cache(self._stage1_warmup)
            self._stage2_cached = self.memory.cache(self._stage2_outliers)
        else:
            self._stage1_cached = self._stage1_warmup
            self._stage2_cached = self._stage2_outliers
        
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
        optuna.logging.set_verbosity(optuna.logging.WARNING)
    
    def run(
        self,
        data: VariantDataset,
        reference: VariantDataset,
        output_dir: Optional[Path] = None,
        generate_plots: bool = True
    ) -> OptimisationResult:
        """
        Run complete three-stage optimisation pipeline.

        Stage 1: Warmup optimisation (if enabled)
        Stage 2: Outlier removal (if enabled)
        Stage 3: Full optimisation

        Parameters
        ----------
        data : VariantDataset
            Input variant data
        reference : VariantDataset
            Reference dataset for regression targets
        output_dir : Path, optional
            Directory to save plots and filtered results
        generate_plots : bool, default=True
            Whether to automatically generate plots after optimization

        Returns
        -------
        OptimisationResult
            Complete optimization results including parameters and scores
        """
        logger.info("="*60)
        logger.info("Starting three-stage optimisation pipeline")
        logger.info("="*60)

        # Calculate regression targets from reference
        logger.info("Calculating regression targets from reference data...")
        targets_by_type = self._calculate_targets(reference)
        logger.info(f"Targets calculated for {len(targets_by_type)} variant types")

        warmup_params = {}
        data_clean = data
        n_removed = 0

        # Stage 1: Warmup
        logger.info("")
        logger.info("="*60)
        logger.info("STAGE 1: WARMUP OPTIMISATION")
        logger.info("="*60)
        if self.config.stage1.enabled:
            logger.info(f"Running warmup with {self.config.stage1.n_trials} trials per variant type...")
            warmup_params = self._run_warmup(data, targets_by_type)
            logger.info(f"✓ Warmup complete for {len(warmup_params)} variant types")
        else:
            logger.info("Stage 1 disabled in configuration")

        # Stage 2: Outlier removal
        logger.info("")
        logger.info("="*60)
        logger.info("STAGE 2: OUTLIER REMOVAL")
        logger.info("="*60)
        if self.config.stage2.enabled and warmup_params:
            logger.info(f"Removing outliers with DNM count range: {self.config.stage2.min_dnm_count}-{self.config.stage2.max_dnm_count}")
            data_clean, n_removed = self._remove_outliers(data, warmup_params)
            logger.info(f"✓ Removed {n_removed} outlier individuals ({len(data_clean)} samples remaining)")
        elif self.config.stage2.enabled and not warmup_params:
            logger.warning("Skipping outlier removal (no warmup parameters available)")
        else:
            logger.info("Stage 2 disabled in configuration")

        # Stage 3: Full optimisation
        logger.info("")
        logger.info("="*60)
        logger.info("STAGE 3: FULL BAYESIAN OPTIMISATION")
        logger.info("="*60)
        if self.config.stage3.enabled:
            logger.info(f"Running full optimisation with {self.config.stage3.n_trials} trials...")
            logger.info(f"Sampler: {self.config.stage3.sampler}, Pruner: {self.config.stage3.pruner}")
            final_params, scores, study = self._run_full_optimisation(data_clean, targets_by_type)
            logger.info(f"✓ Full optimisation complete")
        else:
            logger.warning("Stage 3 disabled in configuration, using warmup parameters as final")
            final_params = warmup_params
            scores = {}
            study = None

        # Create result object
        result = OptimisationResult(
            best_params=final_params,
            best_scores=scores,
            warmup_params=warmup_params if warmup_params else None,
            n_individuals_removed=n_removed,
            study=study
        )

        # Generate plots if requested
        if generate_plots and output_dir and final_params:
            logger.info("Generating optimization result plots")
            try:
                plot_optimization_results(
                    data_df=data_clean.variants,
                    reference_df=reference.variants,
                    best_params=final_params,
                    output_dir=output_dir,
                    save_filtered=True
                )
                save_parameters(final_params, output_dir)
            except Exception as e:
                logger.warning(f"Failed to generate plots: {e}")

        return result
    
    def _calculate_targets(self, reference: VariantDataset) -> Dict[str, np.ndarray]:
        """Calculate regression targets from reference data."""
        targets = {}
        
        for var_type in self.config.optimisation.variant_types:
            ref_subset = reference.filter_by_type(var_type)
            if len(ref_subset) == 0:
                continue
            
            # Count DNMs per sample
            counts = ref_subset.count_by_sample().rename('dnm_count').reset_index()
            
            # Get parental ages
            ages = ref_subset.variants[['SAMPLE', 'paternal_age', 'maternal_age']].drop_duplicates()
            regression_data = ages.merge(counts, on='SAMPLE', how='left').fillna(0)
            
            # Fit regression
            try:
                model = smf.ols(self.config.optimisation.regression_formula, data=regression_data).fit()
                targets[var_type] = model.params.values
                logger.info(f"Calculated targets for {var_type}: {model.params.to_dict()}")
            except Exception as e:
                logger.warning(f"Failed to calculate targets for {var_type}: {e}")
        
        return targets
    
    def _run_warmup(
        self,
        data: VariantDataset,
        targets: Dict[str, np.ndarray]
    ) -> Dict[str, Dict[str, Any]]:
        """Stage 1: Fast warmup optimisation."""
        warmup_params = {}
        
        # Process variant types in parallel if configured
        if self.config.max_workers > 1:
            with ProcessPoolExecutor(max_workers=min(self.config.max_workers, len(targets))) as executor:
                futures = {}
                for var_type in targets:
                    future = executor.submit(
                        self._optimize_variant_type,
                        data.filter_by_type(var_type),
                        targets[var_type],
                        self.config.stage1.n_trials,
                        f"warmup_{var_type}"
                    )
                    futures[var_type] = future

                for var_type, future in futures.items():
                    params, _ = future.result()
                    if params:
                        warmup_params[var_type] = params
        else:
            # Sequential processing for deterministic results
            for var_type in targets:
                params, _ = self._optimize_variant_type(
                    data.filter_by_type(var_type),
                    targets[var_type],
                    self.config.stage1.n_trials,
                    f"warmup_{var_type}"
                )
                if params:
                    warmup_params[var_type] = params
        
        return warmup_params
    
    def _stage1_warmup(self, data_hash: str, config_hash: str) -> Dict[str, Dict[str, Any]]:
        """Cached warmup stage."""
        # This is wrapped by joblib.Memory for caching
        return self._run_warmup(data, targets)
    
    def _remove_outliers(
        self,
        data: VariantDataset,
        warmup_params: Dict[str, Dict[str, Any]]
    ) -> Tuple[VariantDataset, int]:
        """Stage 2: Remove outlier individuals based on filtered DNM counts."""
        
        # Apply warmup filters and count DNMs
        all_filtered = []
        for var_type, params in warmup_params.items():
            var_data = data.filter_by_type(var_type)
            filtered = var_data.apply_filters(params)
            all_filtered.append(filtered.variants)
        
        if not all_filtered:
            return data, 0
        
        # Combine and count total filtered DNMs per individual
        combined = pd.concat(all_filtered, ignore_index=True)
        total_counts = combined.groupby('SAMPLE').size()
        
        logger.info(f"DNM count distribution: min={total_counts.min()}, "
                   f"max={total_counts.max()}, mean={total_counts.mean():.1f}")
        
        # Determine samples to keep
        samples_to_keep = total_counts.index
        
        if self.config.stage2.min_dnm_count:
            samples_to_keep = total_counts[total_counts >= self.config.stage2.min_dnm_count].index
        
        if self.config.stage2.max_dnm_count:
            samples_to_keep = total_counts[
                (total_counts <= self.config.stage2.max_dnm_count) &
                (total_counts.index.isin(samples_to_keep))
            ].index
        
        n_removed = len(total_counts) - len(samples_to_keep)
        
        # Filter data
        filtered_variants = data.variants[data.variants['SAMPLE'].isin(samples_to_keep)].copy()
        
        return VariantDataset(variants=filtered_variants, metadata=data.metadata), n_removed
    
    def _stage2_outliers(self, data_hash: str, warmup_hash: str) -> Tuple[str, int]:
        """Cached outlier removal stage."""
        # Returns hash of cleaned data and number removed
        clean_data, n_removed = self._remove_outliers(data, warmup_params)
        return clean_data.get_hash(), n_removed
    
    def _run_full_optimisation(
        self,
        data: VariantDataset,
        targets: Dict[str, np.ndarray]
    ) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, float], Optional[optuna.Study]]:
        """Stage 3: Full Bayesian optimisation."""
        final_params = {}
        scores = {}
        studies = {}

        for var_type in targets:
            params, best_score = self._optimize_variant_type(
                data.filter_by_type(var_type),
                targets[var_type],
                self.config.stage3.n_trials,
                f"final_{var_type}",
                use_pruning=self.config.stage3.pruner is not None
            )
            if params:
                final_params[var_type] = params
                scores[var_type] = best_score

        return final_params, scores, None  # Return None for study for now
    
    def _optimize_variant_type(
        self,
        data: VariantDataset,
        targets: np.ndarray,
        n_trials: int,
        study_name: str,
        use_pruning: bool = False
    ) -> Tuple[Optional[Dict[str, Any]], float]:
        """Optimize filtering parameters for one variant type.

        Returns:
            Tuple of (best_params, best_score). If optimization fails, returns (None, inf).
        """

        if len(data) == 0:
            return None, float('inf')
        
        # Get parental ages
        parental_info = data.variants[['SAMPLE', 'paternal_age', 'maternal_age']].drop_duplicates()
        
        def objective(trial):
            # Suggest parameters based on data distribution
            params = self._suggest_params(trial, data)
            
            # Apply filters
            filtered = data.apply_filters(params)
            
            if len(filtered) == 0:
                return float('inf')
            
            # Count DNMs
            dnm_counts = filtered.count_by_sample().rename('dnm_count').reset_index()
            regression_data = parental_info.merge(dnm_counts, on='SAMPLE', how='left').fillna(0)

            # Check if we have enough data for regression
            if len(regression_data) < 3:
                logger.debug(f"Not enough samples for regression: {len(regression_data)}")
                return float('inf')

            try:
                # Fit regression
                model = smf.ols(self.config.optimisation.regression_formula, data=regression_data).fit()
                model_params = model.params.values

                # Calculate weighted mean squared error with hybrid metric
                # For small coefficients (|target| < 0.1), use absolute error
                # For large coefficients, use relative error
                # This prevents small coefficients from dominating the loss

                absolute_errors = model_params - targets
                errors = []
                for abs_err, target in zip(absolute_errors, targets):
                    if np.abs(target) < 0.1:
                        # Use absolute error for small coefficients (like insertion slopes)
                        errors.append(abs_err ** 2)
                    else:
                        # Use relative error for large coefficients (like intercepts and SNV slopes)
                        rel_err = abs_err / (np.abs(target) + 1e-10)
                        errors.append(rel_err ** 2)

                errors = np.array(errors)

                # Use intercept_weight for intercept, 1.0 for other coefficients
                intercept_weight = params.get('intercept_weight', 1.0)
                weights = [intercept_weight] + [1.0] * (len(errors) - 1)
                weighted_errors = errors * weights[:len(errors)]
                msre = np.mean(weighted_errors)

                # Report for pruning if enabled
                if use_pruning and trial.number > 0:
                    trial.report(msre, trial.number)
                    if trial.should_prune():
                        raise optuna.TrialPruned()

                return msre

            except Exception as e:
                # Log detailed error information
                logger.warning(f"Regression failed: {e}")
                logger.debug(f"  Regression data shape: {regression_data.shape}")
                logger.debug(f"  DNM count range: {regression_data['dnm_count'].min()}-{regression_data['dnm_count'].max()}")
                logger.debug(f"  Age columns present: {regression_data.columns.tolist()}")
                logger.debug(f"  Formula: {self.config.optimisation.regression_formula}")
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
        
        # Run optimisation
        study.optimize(
            objective,
            n_trials=n_trials,
            n_jobs=1 if self.config.deterministic else self.config.max_workers
        )
        
        logger.info(f"Optimisation for {study_name} finished. Best score: {study.best_value:.6f}")

        # Add detailed logging about the best trial
        try:
            best_params = study.best_params

            # Apply best filters to get final data
            filtered = data.apply_filters(best_params)
            retention_rate = len(filtered) / len(data) * 100 if len(data) > 0 else 0

            # Calculate DNM counts and fit regression
            dnm_counts = filtered.count_by_sample().rename('dnm_count').reset_index()
            regression_data = parental_info.merge(dnm_counts, on='SAMPLE', how='left').fillna(0)

            # Fit regression to get actual coefficients
            model = smf.ols(self.config.optimisation.regression_formula, data=regression_data).fit()
            fitted_coeffs = dict(zip(model.params.index, model.params.values))
            target_coeffs = dict(zip(model.params.index, targets))

            logger.info(f"  Retention: {len(filtered)}/{len(data)} variants ({retention_rate:.1f}%)")
            logger.info(f"  Samples: {len(regression_data)} with DNM counts")
            logger.info(f"  Target coefficients: {', '.join([f'{k}={v:.4f}' for k, v in target_coeffs.items()])}")
            logger.info(f"  Fitted coefficients: {', '.join([f'{k}={v:.4f}' for k, v in fitted_coeffs.items()])}")

            # Calculate and log absolute and relative errors
            errors = {k: fitted_coeffs[k] - target_coeffs[k] for k in target_coeffs.keys()}
            relative_errors = {k: (fitted_coeffs[k] - target_coeffs[k]) / target_coeffs[k] * 100
                             for k in target_coeffs.keys()}
            logger.info(f"  Absolute errors: {', '.join([f'{k}={v:+.4f}' for k, v in errors.items()])}")
            logger.info(f"  Relative errors: {', '.join([f'{k}={v:+.1f}%' for k, v in relative_errors.items()])}")

        except Exception as e:
            logger.debug(f"Could not compute detailed statistics: {e}")

        return study.best_params, study.best_value
    
    def _suggest_params(self, trial: optuna.Trial, data: VariantDataset) -> Dict[str, Any]:
        """Suggest filtering parameters based on data distribution."""
        params = {}

        # Get columns to optimize from configuration
        opt_columns = self.config.optimisation.get_optimisation_columns()
        if not opt_columns:
            return params

        # Check if using CMA-ES sampler (doesn't support discrete parameters)
        is_cmaes = self.config.stage3.sampler == 'cmaes'

        # Build set of columns that are linked (not the first in their group)
        # These should NOT be suggested independently
        linked_groups = self.config.optimisation.get_linked_column_groups()
        skip_columns = set()
        for group in linked_groups:
            if len(group) >= 2:
                # Skip all columns except the first one
                for col_config in group[1:]:
                    skip_columns.add(col_config.name)

        # Suggest parameters for each column using configuration
        for col_config in opt_columns:
            col = col_config.name

            # Skip linked columns (they'll be set to match the first column)
            if col in skip_columns:
                continue

            if col not in data.variants.columns:
                continue

            col_data = pd.to_numeric(data.variants[col], errors='coerce').dropna()
            if len(col_data) == 0:
                continue

            # Use optimization type from configuration
            if col_config.optimisation == 'minimum':
                # Keep LOW values: suggest maximum threshold (filter out values above this)
                lower = col_data.quantile(0.1)
                upper = col_data.max()
                if col_config.dtype == 'int' and not is_cmaes:
                    params[f'max_{col}'] = trial.suggest_int(f'max_{col}', int(lower), int(upper))
                else:
                    # Use float for CMA-ES or if dtype is float
                    value = trial.suggest_float(f'max_{col}', float(lower), float(upper))
                    params[f'max_{col}'] = int(round(value)) if col_config.dtype == 'int' else value

            elif col_config.optimisation == 'maximum':
                # Keep HIGH values: suggest minimum threshold (filter out values below this)
                lower = col_data.min()
                upper = col_data.quantile(0.9)
                if col_config.dtype == 'int' and not is_cmaes:
                    params[f'min_{col}'] = trial.suggest_int(f'min_{col}', int(lower), int(upper))
                else:
                    # Use float for CMA-ES or if dtype is float
                    value = trial.suggest_float(f'min_{col}', float(lower), float(upper))
                    params[f'min_{col}'] = int(round(value)) if col_config.dtype == 'int' else value

            elif col_config.optimisation == 'range':
                # Range: handle symmetric vs regular range constraints
                if col_config.range_constraint:
                    # Check if it's a symmetric range constraint (has 'lower' and 'scale')
                    if hasattr(col_config.range_constraint, 'lower') and hasattr(col_config.range_constraint, 'scale'):
                        # Symmetric range: suggest only lower, calculate upper as scale - lower
                        lower_bound = col_config.range_constraint.lower
                        upper_bound = col_config.range_constraint.scale - col_config.range_constraint.lower

                        if col_config.dtype == 'int' and not is_cmaes:
                            suggested_lower = trial.suggest_int(f'min_{col}', int(lower_bound), int(upper_bound))
                        else:
                            value = trial.suggest_float(f'min_{col}', float(lower_bound), float(upper_bound))
                            suggested_lower = int(round(value)) if col_config.dtype == 'int' else value

                        # Calculate symmetric upper bound
                        suggested_upper = col_config.range_constraint.scale - suggested_lower
                        params[f'min_{col}'] = suggested_lower
                        params[f'max_{col}'] = suggested_upper
                    else:
                        # Regular range: suggest both min and max independently
                        lower_bound = col_config.range_constraint.min
                        upper_bound = col_config.range_constraint.max

                        if col_config.dtype == 'int' and not is_cmaes:
                            params[f'min_{col}'] = trial.suggest_int(f'min_{col}', int(lower_bound), int(upper_bound))
                            params[f'max_{col}'] = trial.suggest_int(f'max_{col}', int(lower_bound), int(upper_bound))
                        else:
                            value_min = trial.suggest_float(f'min_{col}', float(lower_bound), float(upper_bound))
                            value_max = trial.suggest_float(f'max_{col}', float(lower_bound), float(upper_bound))
                            if col_config.dtype == 'int':
                                params[f'min_{col}'] = int(round(value_min))
                                params[f'max_{col}'] = int(round(value_max))
                            else:
                                params[f'min_{col}'] = value_min
                                params[f'max_{col}'] = value_max
                else:
                    # No constraint: use data bounds
                    lower_bound = col_data.min()
                    upper_bound = col_data.max()

                    if col_config.dtype == 'int' and not is_cmaes:
                        params[f'min_{col}'] = trial.suggest_int(f'min_{col}', int(lower_bound), int(upper_bound))
                        params[f'max_{col}'] = trial.suggest_int(f'max_{col}', int(lower_bound), int(upper_bound))
                    else:
                        value_min = trial.suggest_float(f'min_{col}', float(lower_bound), float(upper_bound))
                        value_max = trial.suggest_float(f'max_{col}', float(lower_bound), float(upper_bound))
                        if col_config.dtype == 'int':
                            params[f'min_{col}'] = int(round(value_min))
                            params[f'max_{col}'] = int(round(value_max))
                        else:
                            params[f'min_{col}'] = value_min
                            params[f'max_{col}'] = value_max

        # Handle linked columns - copy first column's value to linked columns
        for group in linked_groups:
            if len(group) < 2:
                continue

            # Find parameters for the first column in the group
            first_col = group[0].name
            param_keys = [k for k in params.keys() if k.endswith(f'_{first_col}')]

            if param_keys:
                # Copy the first column's value to all linked columns
                for key in param_keys:
                    param_value = params[key]
                    prefix = key.split('_')[0]  # 'min' or 'max'

                    # Apply same value to all linked columns
                    for linked_col in group[1:]:
                        linked_key = f'{prefix}_{linked_col.name}'
                        params[linked_key] = param_value

        # Add tunable intercept weight for regression objective
        # Higher weights (0.5-2.0) give intercept more importance relative to slope
        params['intercept_weight'] = trial.suggest_float('intercept_weight', 0.5, 2.0)

        return params
