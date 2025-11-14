"""Command-line interface for variant optimizer."""

import typer
from pathlib import Path
from typing import Optional, List
import yaml
import logging
import sys

from .api import optimize_filters, run_optimization, StageRunner
from .config import PipelineConfig, PRESETS
from .data import VariantDataset


app = typer.Typer(
    name="variant-optimize",
    help="Three-stage Bayesian optimization for genomic variant filtering",
    add_completion=False
)


def setup_logging(verbose: int):
    """Configure logging based on verbosity level."""
    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    level = levels[min(verbose, 2)]
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )


@app.command()
def run(
    data: Path = typer.Argument(..., help="Path to variant data TSV file"),
    reference: Path = typer.Argument(..., help="Path to reference data TSV file"),
    config: Optional[Path] = typer.Option(None, "--config", "-c", help="Configuration file (YAML)"),
    preset: Optional[str] = typer.Option(
        "balanced", "--preset", "-p",
        help=f"Configuration preset: {', '.join(PRESETS.keys())}"
    ),
    output: Path = typer.Option("./results", "--output", "-o", help="Output directory"),
    n_trials: Optional[int] = typer.Option(None, "--n-trials", "-n", help="Number of optimization trials"),
    min_dnm: Optional[int] = typer.Option(None, "--min-dnm", help="Minimum DNM count for outlier removal"),
    max_dnm: Optional[int] = typer.Option(None, "--max-dnm", help="Maximum DNM count for outlier removal"),
    workers: Optional[int] = typer.Option(None, "--workers", "-w", help="Number of parallel workers"),
    seed: int = typer.Option(42, "--seed", "-s", help="Random seed for reproducibility"),
    override: Optional[List[str]] = typer.Option(
        None, "--override", "-O",
        help="Configuration overrides (e.g., -O stage1.n_trials=100 -O stage3.sampler=cmaes)"
    ),
    verbose: int = typer.Option(0, "--verbose", "-v", count=True, help="Increase verbosity"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show configuration without running")
):
    """
    Run three-stage variant filtering optimization.
    
    Examples:
    
        # Basic usage with default settings
        variant-optimize data.tsv reference.tsv
        
        # Fast testing
        variant-optimize data.tsv reference.tsv --preset fast
        
        # Thorough optimization with custom trials
        variant-optimize data.tsv reference.tsv --preset thorough --n-trials 1000
        
        # Custom outlier removal thresholds
        variant-optimize data.tsv reference.tsv --min-dnm 10 --max-dnm 200
        
        # Use configuration file with overrides
        variant-optimize data.tsv reference.tsv -c config.yaml -O stage3.n_trials=1000
        
        # Parallel processing with 8 workers
        variant-optimize data.tsv reference.tsv --workers 8
    """
    setup_logging(verbose)
    
    # Validate inputs
    if not data.exists():
        typer.echo(f"Error: Data file {data} not found", err=True)
        raise typer.Exit(1)
    
    if not reference.exists():
        typer.echo(f"Error: Reference file {reference} not found", err=True)
        raise typer.Exit(1)
    
    # Load or create configuration
    if config and config.exists():
        pipeline_config = PipelineConfig.from_yaml(config)
        typer.echo(f"Loaded configuration from {config}")
    elif preset:
        pipeline_config = PipelineConfig.from_preset(preset)
        typer.echo(f"Using preset: {preset}")
    else:
        pipeline_config = PipelineConfig()
        typer.echo("Using default configuration")
    
    # Apply command-line overrides
    if n_trials:
        pipeline_config.stage3.n_trials = n_trials
    if min_dnm is not None:
        pipeline_config.stage2.min_dnm_count = min_dnm
    if max_dnm is not None:
        pipeline_config.stage2.max_dnm_count = max_dnm
    if workers:
        pipeline_config.max_workers = workers
    
    pipeline_config.seed = seed
    
    # Apply string overrides
    if override:
        for override_str in override:
            if '=' not in override_str:
                typer.echo(f"Error: Invalid override format: {override_str}", err=True)
                typer.echo("Expected format: key=value (e.g., stage1.n_trials=100)", err=True)
                raise typer.Exit(1)
            
            key, value = override_str.split('=', 1)
            keys = key.split('.')
            
            # Navigate to the nested field
            obj = pipeline_config
            for k in keys[:-1]:
                obj = getattr(obj, k)
            
            # Set the value (attempt type conversion)
            try:
                # Try to evaluate as Python literal
                import ast
                parsed_value = ast.literal_eval(value)
            except (ValueError, SyntaxError):
                # Keep as string
                parsed_value = value
            
            setattr(obj, keys[-1], parsed_value)
    
    # Show configuration if verbose or dry run
    if verbose or dry_run:
        typer.echo("\nConfiguration:")
        typer.echo("-" * 40)
        config_dict = pipeline_config.model_dump()
        typer.echo(yaml.dump(config_dict, default_flow_style=False))
        typer.echo("-" * 40)
    
    if dry_run:
        typer.echo("\nDry run - not executing optimization")
        raise typer.Exit(0)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Save configuration
    pipeline_config.to_yaml(output / "config.yaml")
    
    # Run optimization
    typer.echo(f"\nStarting optimization with {pipeline_config.stage3.n_trials} trials...")
    typer.echo(f"Output directory: {output}")
    
    try:
        result = run_optimization(
            data=data,
            reference=reference,
            config=pipeline_config
        )
        
        # Save results
        with open(output / "optimal_params.yaml", 'w') as f:
            yaml.dump(result.best_params, f, default_flow_style=False)
        
        with open(output / "summary.txt", 'w') as f:
            f.write(result.summary)
        
        # Display summary
        typer.echo("\n" + "="*50)
        typer.echo("OPTIMIZATION COMPLETE")
        typer.echo("="*50)
        typer.echo(result.summary)
        typer.echo(f"\nResults saved to: {output}")
        
    except Exception as e:
        typer.echo(f"\nError during optimization: {e}", err=True)
        if verbose > 0:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)


@app.command()
def show_presets():
    """Show available configuration presets."""
    typer.echo("Available presets:\n")
    
    for name, config in PRESETS.items():
        typer.echo(f"{name}:")
        typer.echo(f"  Stage 1: {config.stage1.n_trials} warmup trials")
        typer.echo(f"  Stage 2: outlier removal {'enabled' if config.stage2.enabled else 'disabled'}")
        if config.stage2.enabled and config.stage2.min_dnm_count:
            typer.echo(f"    - Min DNMs: {config.stage2.min_dnm_count}")
            typer.echo(f"    - Max DNMs: {config.stage2.max_dnm_count}")
        typer.echo(f"  Stage 3: {config.stage3.n_trials} optimization trials")
        typer.echo(f"    - Sampler: {config.stage3.sampler}")
        typer.echo(f"  Workers: {config.max_workers}")
        typer.echo()


@app.command()
def validate(
    data: Path = typer.Argument(..., help="Path to variant data TSV file"),
    reference: Optional[Path] = typer.Option(None, help="Path to reference data TSV file"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show detailed information")
):
    """
    Validate data files and show statistics.
    
    Useful for checking data before running optimization.
    """
    typer.echo(f"Validating {data}...")
    
    try:
        dataset = VariantDataset.from_tsv(data)
        summary = dataset.summary
        
        typer.echo(f"✓ Successfully loaded {summary['n_variants']} variants")
        typer.echo(f"  Samples: {summary['n_samples']}")
        typer.echo(f"  Chromosomes: {summary['chromosomes']}")
        
        if summary['variant_types']:
            typer.echo("  Variant types:")
            for vtype, count in summary['variant_types'].items():
                typer.echo(f"    - {vtype}: {count}")
        
        if verbose:
            typer.echo(f"\n  Columns: {', '.join(dataset.variants.columns)}")
            
            # Show data types
            typer.echo("\n  Data types:")
            for col, dtype in dataset.variants.dtypes.items():
                typer.echo(f"    - {col}: {dtype}")
            
            # Memory usage
            memory_mb = dataset.variants.memory_usage(deep=True).sum() / 1024**2
            typer.echo(f"\n  Memory usage: {memory_mb:.2f} MB")
        
        if reference:
            typer.echo(f"\nValidating reference {reference}...")
            ref_dataset = VariantDataset.from_tsv(reference)
            ref_summary = ref_dataset.summary
            
            typer.echo(f"✓ Successfully loaded {ref_summary['n_variants']} reference variants")
            typer.echo(f"  Samples: {ref_summary['n_samples']}")
            
            # Check compatibility
            data_types = set(summary['variant_types'].keys())
            ref_types = set(ref_summary['variant_types'].keys())
            common_types = data_types & ref_types
            
            if common_types:
                typer.echo(f"\n✓ Compatible variant types: {', '.join(common_types)}")
            else:
                typer.echo("\n⚠ Warning: No common variant types between data and reference", err=True)
        
        typer.echo("\n✓ Validation complete - files are ready for optimization")
        
    except Exception as e:
        typer.echo(f"\n✗ Validation failed: {e}", err=True)
        raise typer.Exit(1)


@app.command()
def init(
    output: Path = typer.Option(Path("config.yaml"), "--output", "-o", help="Output configuration file"),
    preset: str = typer.Option("balanced", "--preset", "-p", help="Base preset to use"),
    interactive: bool = typer.Option(False, "--interactive", "-i", help="Interactive configuration")
):
    """
    Initialize a configuration file.
    
    Creates a YAML configuration file that can be edited and used with --config.
    """
    config = PipelineConfig.from_preset(preset)
    
    if interactive:
        typer.echo("Interactive configuration setup")
        typer.echo("-" * 40)
        
        # Stage 1
        if typer.confirm("Configure warmup stage?", default=True):
            config.stage1.enabled = typer.confirm("  Enable warmup?", default=True)
            if config.stage1.enabled:
                config.stage1.n_trials = typer.prompt("  Warmup trials", default=50, type=int)
        
        # Stage 2
        if typer.confirm("Configure outlier removal?", default=True):
            config.stage2.enabled = typer.confirm("  Enable outlier removal?", default=True)
            if config.stage2.enabled:
                config.stage2.min_dnm_count = typer.prompt("  Minimum DNM count", default=5, type=int)
                config.stage2.max_dnm_count = typer.prompt("  Maximum DNM count", default=300, type=int)
        
        # Stage 3
        if typer.confirm("Configure optimization stage?", default=True):
            config.stage3.n_trials = typer.prompt("  Optimization trials", default=500, type=int)
            sampler = typer.prompt(
                "  Sampler (tpe/cmaes/random)",
                default="tpe",
                type=typer.Choice(["tpe", "cmaes", "random"])
            )
            config.stage3.sampler = sampler
        
        # Performance
        config.max_workers = typer.prompt("Number of parallel workers", default=4, type=int)
        config.seed = typer.prompt("Random seed", default=42, type=int)
    
    # Save configuration
    config.to_yaml(output)
    typer.echo(f"\nConfiguration saved to {output}")
    typer.echo(f"Use with: variant-optimize data.tsv reference.tsv --config {output}")


if __name__ == "__main__":
    app()
