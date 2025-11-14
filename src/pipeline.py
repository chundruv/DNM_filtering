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

from .config import PipelineConfig
from .data import VariantDataset


logger = logging.getLogger(__name__)


@dataclass
class OptimizationResult:
    """Container for optimization results."""
    best_params: Dict[str, Dict[str, Any]]
    best_scores: Dict[str, float]
    warmup_params: Optional[Dict[str, Dict[str, Any]]]
    n_individuals_removed: int
    study: Optional[optuna.Study]
    
    @property
    def summary(self) -> str:
        """Human-readable summary."""
        lines = ["Optimization Results:"]
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


class OptimizationPipeline:
    """Three-stage optimization pipeline with caching and parallelization."""
    
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
        reference: VariantDataset
    ) -> OptimizationResult:
        """
        Run complete three-stage optimization pipeline.
        
        Stage 1: Warmup optimization (if enabled)
        Stage 2: Outlier removal (if enabled)
        Stage 3: Full optimization
        """
        logger.info("Starting three-stage optimization pipeline")
        
        # Calculate regression targets from reference
        targets_by_type = self._calculate_targets(reference)
        
        warmup_params = {}
        data_clean = data
        n_removed = 0
        
        # Stage 1: Warmup
        if self.config.stage1.enabled:
            logger.info("Stage 1: Running warmup optimization")
            warmup_params = self._run_warmup(data, targets_by_type)
            logger.info(f"Warmup complete for {len(warmup_params)} variant types")
        
        # Stage 2: Outlier removal
        if self.config.stage2.enabled and warmup_params:
            logger.info("Stage 2: Removing outliers based on filtered DNM counts")
            data_clean, n_removed = self._remove_outliers(data, warmup_params)
            logger.info(f"Removed {n_removed} outlier individuals")
        
        # Stage 3: Full optimization
        logger.info("Stage 3: Running full optimization")
        final_params, scores, study = self._run_full_optimization(data_clean, targets_by_type)
        
        return OptimizationResult(
            best_params=final_params,
            best_scores=scores,
            warmup_params=warmup_params if warmup_params else None,
            n_individuals_removed=n_removed,
            study=study
        )
    
    def _calculate_targets(self, reference: VariantDataset) -> Dict[str, np.ndarray]:
        """Calculate regression targets from reference data."""
        targets = {}
        
        for var_type in self.config.optimization.variant_types:
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
                model = smf.ols(self.config.optimization.regression_formula, data=regression_data).fit()
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
        """Stage 1: Fast warmup optimization."""
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
                    result = future.result()
                    if result:
                        warmup_params[var_type] = result
        else:
            # Sequential processing for deterministic results
            for var_type in targets:
                result = self._optimize_variant_type(
                    data.filter_by_type(var_type),
                    targets[var_type],
                    self.config.stage1.n_trials,
                    f"warmup_{var_type}"
                )
                if result:
                    warmup_params[var_type] = result
        
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
    
    def _run_full_optimization(
        self,
        data: VariantDataset,
        targets: Dict[str, np.ndarray]
    ) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, float], Optional[optuna.Study]]:
        """Stage 3: Full Bayesian optimization."""
        final_params = {}
        scores = {}
        studies = {}
        
        for var_type in targets:
            params = self._optimize_variant_type(
                data.filter_by_type(var_type),
                targets[var_type],
                self.config.stage3.n_trials,
                f"final_{var_type}",
                use_pruning=self.config.stage3.pruner is not None
            )
            if params:
                final_params[var_type] = params
                # Store score (we'd need to modify _optimize_variant_type to return this)
                scores[var_type] = 0.0  # Placeholder
        
        return final_params, scores, None  # Return None for study for now
    
    def _optimize_variant_type(
        self,
        data: VariantDataset,
        targets: np.ndarray,
        n_trials: int,
        study_name: str,
        use_pruning: bool = False
    ) -> Optional[Dict[str, Any]]:
        """Optimize filtering parameters for one variant type."""
        
        if len(data) == 0:
            return None
        
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
            
            try:
                # Fit regression
                model = smf.ols(self.config.optimization.regression_formula, data=regression_data).fit()
                model_params = model.params.values
                
                # Calculate weighted MSE
                squared_errors = (model_params - targets) ** 2
                weights = self.config.optimization.regression_weights
                if len(weights) < len(squared_errors):
                    weights = weights + [1.0] * (len(squared_errors) - len(weights))
                weighted_errors = squared_errors * weights[:len(squared_errors)]
                mse = np.mean(weighted_errors)
                
                # Report for pruning if enabled
                if use_pruning and trial.number > 0:
                    trial.report(mse, trial.number)
                    if trial.should_prune():
                        raise optuna.TrialPruned()
                
                return mse
                
            except Exception as e:
                logger.debug(f"Regression failed: {e}")
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
        
        # Run optimization
        study.optimize(
            objective,
            n_trials=n_trials,
            n_jobs=1 if self.config.deterministic else self.config.max_workers
        )
        
        logger.info(f"Optimization for {study_name} finished. Best score: {study.best_value:.6f}")
        
        return study.best_params
    
    def _suggest_params(self, trial: optuna.Trial, data: VariantDataset) -> Dict[str, Any]:
        """Suggest filtering parameters based on data distribution."""
        params = {}
        
        # Get columns to optimize
        cols_to_optimize = self.config.optimization.column_names
        if cols_to_optimize is None:
            # Use defaults based on available columns
            possible_cols = ['VAF', 'child_DP', 'father_DP', 'mother_DP', 'QUAL', 'GQ']
            cols_to_optimize = [c for c in possible_cols if c in data.variants.columns]
        
        # Suggest parameters for each column
        for col in cols_to_optimize:
            if col not in data.variants.columns:
                continue
            
            col_data = pd.to_numeric(data.variants[col], errors='coerce').dropna()
            if len(col_data) == 0:
                continue
            
            # Determine parameter type and bounds
            if 'DP' in col or 'depth' in col.lower() or 'coverage' in col.lower():
                # Depth/coverage: suggest minimum threshold
                lower = col_data.min()
                upper = col_data.quantile(0.75)
                if col_data.dtype.kind == 'i':  # Integer
                    params[f'min_{col}'] = trial.suggest_int(f'min_{col}', int(lower), int(upper))
                else:
                    params[f'min_{col}'] = trial.suggest_float(f'min_{col}', lower, upper)
            
            elif 'VAF' in col or 'freq' in col.lower():
                # Allele frequency: suggest symmetric range
                if col_data.max() <= 1.0:
                    # 0-1 scale
                    params[f'min_{col}'] = trial.suggest_float(f'min_{col}', 0.15, 0.4)
                else:
                    # 0-100 scale
                    params[f'min_{col}'] = trial.suggest_float(f'min_{col}', 15, 40)
            
            elif 'QUAL' in col or 'quality' in col.lower() or 'GQ' in col:
                # Quality scores: minimum threshold
                lower = col_data.quantile(0.1)
                upper = col_data.quantile(0.9)
                params[f'min_{col}'] = trial.suggest_float(f'min_{col}', lower, upper)
        
        # Link parent depths if both present
        if 'min_father_DP' in params and 'min_mother_DP' in params:
            # Use same threshold for both parents
            parent_threshold = params['min_father_DP']
            params['min_mother_DP'] = parent_threshold
        
        return params
