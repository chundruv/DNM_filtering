"""Progressive disclosure API for variant optimisation."""

from pathlib import Path
from typing import Dict, Any, Optional, Union
import logging

from .config import PipelineConfig, load_config
from .data import VariantDataset
from .pipeline import OptimisationPipeline, MemoryEfficientPipeline, OptimisationResult


logger = logging.getLogger(__name__)


# Level 1: Simple function interface (80% of use cases)
def optimize_filters(
    data_path: Union[str, Path],
    reference_path: Union[str, Path],
    n_trials: int = 100,
    preset: str = "balanced",
    output_dir: Optional[Union[str, Path]] = None,
    seed: int = 42,
    memory_efficient: bool = False
) -> Dict[str, Dict[str, Any]]:
    """
    Simple interface for variant filtering optimisation.
    
    This is the easiest way to use the optimizer - just provide your data
    and reference files, and it handles everything else.
    
    Parameters
    ----------
    data_path : str or Path
        Path to variant data (TSV or Parquet)
    reference_path : str or Path
        Path to reference data (TSV or Parquet)
    n_trials : int, default=100
        Number of Stage 3 trials (overrides the preset's value)
    preset : str, default="balanced"
        Configuration preset: "fast", "balanced", "thorough", or "noisy_data"
    output_dir : str or Path, optional
        Directory to save results
    seed : int, default=42
        Random seed for reproducibility
    memory_efficient : bool, default=False
        Currently has no effect (OptimisationPipeline and MemoryEfficientPipeline are the same class)
    
    Returns
    -------
    dict
        Optimized filtering parameters for each variant type
    
    Examples
    --------
    >>> # Basic usage
    >>> params = optimize_filters("data.tsv", "reference.tsv")
    
    >>> # Fast testing
    >>> params = optimize_filters("data.tsv", "reference.tsv", preset="fast")
    
    >>> # Thorough optimisation
    >>> params = optimize_filters("data.tsv", "reference.tsv", n_trials=500, preset="thorough")
    
    >>> # Memory-efficient mode for large datasets
    >>> params = optimize_filters("data.tsv", "reference.tsv", memory_efficient=True)
    """
    # Load preset configuration
    config = PipelineConfig.from_preset(preset)
    config.stage3.n_trials = n_trials
    config.seed = seed
    
    # Load data
    logger.info(f"Loading data from {data_path}")
    data = VariantDataset.from_path(
        Path(data_path),
        sample_col=config.optimisation.sample_id_column_data,
        paternal_age_col=config.optimisation.paternal_age_column_data,
        maternal_age_col=config.optimisation.maternal_age_column_data,
        reference_col=config.optimisation.reference_column_data,
        alternate_col=config.optimisation.alternate_column_data,
        chrom_col=config.optimisation.chromosome_column_data,
        pos_col=config.optimisation.position_column_data,
        columns=config.optimisation.get_input_columns('data'),
    )

    logger.info(f"Loading reference from {reference_path}")
    reference = VariantDataset.from_path(
        Path(reference_path),
        sample_col=config.optimisation.sample_id_column_reference,
        paternal_age_col=config.optimisation.paternal_age_column_reference,
        maternal_age_col=config.optimisation.maternal_age_column_reference,
        reference_col=config.optimisation.reference_column_reference,
        alternate_col=config.optimisation.alternate_column_reference,
        chrom_col=config.optimisation.chromosome_column_reference,
        pos_col=config.optimisation.position_column_reference,
        columns=config.optimisation.get_input_columns('reference'),
    )
    
    # Run optimisation (use memory-efficient pipeline if requested)
    if memory_efficient:
        pipeline = MemoryEfficientPipeline(config)
    else:
        pipeline = OptimisationPipeline(config)

    # Convert output_dir to Path if provided
    output_path = Path(output_dir) if output_dir else None

    # Run pipeline with automatic plotting
    result = pipeline.run(data, reference, output_dir=output_path, generate_plots=True)

    # Save additional results if requested
    if output_path:
        output_path.mkdir(parents=True, exist_ok=True)

        # Save parameters as YAML
        import yaml
        with open(output_path / "optimal_params.yaml", 'w') as f:
            yaml.dump(result.best_params, f)

        # Save summary
        with open(output_path / "summary.txt", 'w') as f:
            f.write(result.summary)

        logger.info(f"Results saved to {output_path}")

    return result.best_params


# Level 2: Configuration-based interface (advanced users)
def run_optimisation(
    data: Union[VariantDataset, Path, str],
    reference: Union[VariantDataset, Path, str],
    config: Optional[PipelineConfig] = None,
    config_file: Optional[Path] = None,
    output_dir: Optional[Union[str, Path]] = None,
    generate_plots: bool = True,
    memory_efficient: bool = False,
    **kwargs
) -> OptimisationResult:
    """
    Run optimisation with detailed configuration control.

    This interface provides more control over the optimisation process
    while still handling the pipeline orchestration.

    Parameters
    ----------
    data : VariantDataset, Path, or str
        Variant data (dataset object or path to file)
    reference : VariantDataset, Path, or str
        Reference data (dataset object or path to file)
    config : PipelineConfig, optional
        Configuration object. If None, built from preset/config_file/kwargs
    config_file : Path, optional
        Path to YAML configuration file
    output_dir : str or Path, optional
        Directory to save plots and filtered results
    generate_plots : bool, default=True
        Whether to automatically generate plots after optimization
    memory_efficient : bool, default=False
        Currently has no effect (OptimisationPipeline and MemoryEfficientPipeline are the same class)
    **kwargs
        `preset` plus configuration overrides (stage3__n_trials=1000 -> stage3.n_trials).
        Ignored when `config` is given.

    Returns
    -------
    OptimisationResult
        Complete optimisation results including parameters, scores, and metadata

    Examples
    --------
    >>> # With custom configuration
    >>> config = PipelineConfig(
    ...     stage1=Stage1Config(n_trials=100),
    ...     stage2=Stage2Config(min_dnm_count=10, max_dnm_count=200),
    ...     stage3=Stage3Config(n_trials=1000, sampler="cmaes")
    ... )
    >>> result = run_optimisation("data.tsv", "reference.tsv", config=config)

    >>> # From configuration file
    >>> result = run_optimisation("data.tsv", "reference.tsv", config_file="custom.yaml")

    >>> # With overrides and automatic plotting
    >>> result = run_optimisation(
    ...     "data.tsv", "reference.tsv",
    ...     preset="balanced",
    ...     max_workers=8,
    ...     stage3__n_trials=1000,
    ...     output_dir="results/"
    ... )
    
    >>> # Memory-efficient mode for large datasets
    >>> result = run_optimisation("data.tsv", "reference.tsv", memory_efficient=True)
    """
    # Load configuration
    if config is None:
        # Process kwargs for configuration overrides
        preset = kwargs.pop('preset', None)
        cli_overrides = {}

        # Convert kwargs to dot notation for nested config
        for key, value in kwargs.items():
            # Convert __ to . for nested access (e.g., stage3__n_trials -> stage3.n_trials)
            key = key.replace('__', '.')
            cli_overrides[key] = value

        config = load_config(
            config_file=config_file,
            preset=preset,
            cli_overrides=cli_overrides if cli_overrides else None
        )

    # Load data if needed
    if isinstance(data, (str, Path)):
        data = VariantDataset.from_path(
            Path(data),
            sample_col=config.optimisation.sample_id_column_data,
            paternal_age_col=config.optimisation.paternal_age_column_data,
            maternal_age_col=config.optimisation.maternal_age_column_data,
            reference_col=config.optimisation.reference_column_data,
            alternate_col=config.optimisation.alternate_column_data,
            chrom_col=config.optimisation.chromosome_column_data,
            pos_col=config.optimisation.position_column_data,
            columns=config.optimisation.get_input_columns('data'),
        )

    if isinstance(reference, (str, Path)):
        reference = VariantDataset.from_path(
            Path(reference),
            sample_col=config.optimisation.sample_id_column_reference,
            paternal_age_col=config.optimisation.paternal_age_column_reference,
            maternal_age_col=config.optimisation.maternal_age_column_reference,
            reference_col=config.optimisation.reference_column_reference,
            alternate_col=config.optimisation.alternate_column_reference,
            chrom_col=config.optimisation.chromosome_column_reference,
            pos_col=config.optimisation.position_column_reference,
            columns=config.optimisation.get_input_columns('reference'),
        )

    # Run pipeline with automatic plotting (use memory-efficient pipeline if requested)
    if memory_efficient:
        pipeline = MemoryEfficientPipeline(config)
    else:
        pipeline = OptimisationPipeline(config)
    output_path = Path(output_dir) if output_dir else None
    return pipeline.run(data, reference, output_dir=output_path, generate_plots=generate_plots)


# Level 3: Direct stage control (researchers, method development)
class StageRunner:
    """
    Direct control over individual pipeline stages.

    Runs the same stages as `run_optimisation`, one at a time, so intermediate
    results can be inspected or modified between stages.

    Examples
    --------
    >>> runner = StageRunner(PipelineConfig.from_yaml("UKB_filter.yaml"))
    >>> targets = runner.calculate_targets(reference)
    >>> warmup_params = runner.warmup(data, targets, n_trials=100)
    >>> clean_data = runner.remove_outliers(data, warmup_params, min_dnm=5, max_dnm=200)
    >>> final_params = runner.optimize(clean_data, targets, n_trials=1000)
    """

    def __init__(self, config: Optional[PipelineConfig] = None):
        """Initialize with optional base configuration."""
        self.config = config or PipelineConfig()
        self.pipeline = OptimisationPipeline(self.config)

    def calculate_targets(self, reference: VariantDataset) -> Dict[str, Any]:
        """Fit dnm_count ~ paternal_age on the reference, per variant type."""
        return self.pipeline._calculate_targets(reference)

    def warmup(
        self,
        data: VariantDataset,
        targets: Dict[str, Any],
        n_trials: int = 50,
    ) -> Dict[str, Dict[str, Any]]:
        """Run the warmup stage (Stage 1) with `n_trials` trials per variant type."""
        original_trials = self.config.stage1.n_trials
        self.config.stage1.n_trials = n_trials

        try:
            result = self.pipeline._run_warmup_sequential(data, targets)
        finally:
            self.config.stage1.n_trials = original_trials

        return result

    def remove_outliers(
        self,
        data: VariantDataset,
        warmup_params: Dict[str, Dict[str, Any]],
        min_dnm: Optional[int] = None,
        max_dnm: Optional[int] = None,
    ) -> VariantDataset:
        """Run outlier removal (Stage 2) and return a dataset of the kept samples.

        min_dnm / max_dnm override the configured thresholds for this call only.
        """
        original = (self.config.stage2.min_dnm_count, self.config.stage2.max_dnm_count)
        if min_dnm is not None:
            self.config.stage2.min_dnm_count = min_dnm
        if max_dnm is not None:
            self.config.stage2.max_dnm_count = max_dnm

        try:
            samples_to_keep, n_removed = self.pipeline._identify_outliers(data, warmup_params)
        finally:
            self.config.stage2.min_dnm_count, self.config.stage2.max_dnm_count = original

        logger.info(f"Removed {n_removed} outliers")
        kept = data.variants[data.variants['SAMPLE'].isin(samples_to_keep)]
        return VariantDataset(variants=kept.copy(), metadata=dict(data.metadata))

    def optimize(
        self,
        data: VariantDataset,
        targets: Dict[str, Any],
        n_trials: int = 500,
    ) -> Dict[str, Dict[str, Any]]:
        """Run full optimisation (Stage 3), once, with `n_trials` trials per variant type."""
        original_trials = self.config.stage3.n_trials
        self.config.stage3.n_trials = n_trials

        try:
            params, _ = self.pipeline._run_full_optimisation_sequential(data, targets)
        finally:
            self.config.stage3.n_trials = original_trials

        return params


# Convenience function for calculating targets
def calculate_regression_targets(
    reference: Union[VariantDataset, Path, str],
    variant_types: Optional[list] = None,
    formula: str = "dnm_count ~ paternal_age + maternal_age"
) -> Dict[str, Any]:
    """
    Calculate regression targets from reference data.
    
    Useful for custom pipelines or debugging.
    
    Parameters
    ----------
    reference : VariantDataset or path
        Reference dataset
    variant_types : list, optional
        Variant types to process. Default: ["SNV", "Insertion", "Deletion"]
    formula : str
        Currently ignored; targets are always fitted with dnm_count ~ paternal_age.
        Reference columns use the default names (sample_id, paternal_age, ref, alt).
    
    Returns
    -------
    dict
        Regression coefficients for each variant type
    """
    if variant_types is None:
        variant_types = ["SNV", "Insertion", "Deletion"]

    config = PipelineConfig()

    if isinstance(reference, (str, Path)):
        reference = VariantDataset.from_path(
            Path(reference),
            sample_col=config.optimisation.sample_id_column_reference,
            paternal_age_col=config.optimisation.paternal_age_column_reference,
            maternal_age_col=config.optimisation.maternal_age_column_reference,
            reference_col=config.optimisation.reference_column_reference,
            alternate_col=config.optimisation.alternate_column_reference,
            chrom_col=config.optimisation.chromosome_column_reference,
            pos_col=config.optimisation.position_column_reference,
            columns=config.optimisation.get_input_columns('reference'),
        )
    config.optimisation.variant_types = variant_types
    config.optimisation.regression_formula = formula
    
    pipeline = OptimisationPipeline(config)
    return pipeline._calculate_targets(reference)
