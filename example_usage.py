#!/usr/bin/env python3
"""
Example: Using the refactored three-stage optimization pipeline.

This demonstrates the improvements:
- Progressive API (3 levels)
- Performance optimizations
- Configuration management
- Reproducibility
"""

import sys
sys.path.insert(0, 'src')

from variant_optimizer import (
    optimize_filters,
    run_optimization,
    PipelineConfig,
    VariantDataset,
    StageRunner
)
from pathlib import Path

# Create sample data files for testing (normally you'd have real data)
def create_sample_data():
    """Create sample TSV files for demonstration."""
    import pandas as pd
    import numpy as np
    
    np.random.seed(42)
    n_samples = 100
    n_variants = 1000
    
    # Create sample variant data
    data = pd.DataFrame({
        'SAMPLE': np.repeat(range(n_samples), n_variants // n_samples),
        'CHROM': np.random.choice(['chr1', 'chr2', 'chr3'], n_variants),
        'POS': np.random.randint(1000, 100000, n_variants),
        'REF': np.random.choice(['A', 'C', 'G', 'T'], n_variants),
        'ALT': np.random.choice(['A', 'C', 'G', 'T'], n_variants),
        'QUAL': np.random.uniform(10, 60, n_variants),
        'DP': np.random.randint(5, 100, n_variants),
        'VAF': np.random.uniform(0.1, 0.9, n_variants),
        'child_DP': np.random.randint(10, 80, n_variants),
        'father_DP': np.random.randint(10, 80, n_variants),
        'mother_DP': np.random.randint(10, 80, n_variants),
        'GQ': np.random.uniform(20, 99, n_variants),
        'paternal_age': np.random.uniform(25, 45, n_variants),
        'maternal_age': np.random.uniform(23, 42, n_variants)
    })
    
    # Save data
    data.to_csv('sample_data.tsv', sep='\t', index=False)
    
    # Create reference with similar structure
    reference = data.sample(n=500).copy()
    reference['SAMPLE'] = np.random.choice(range(50), 500)
    reference.to_csv('sample_reference.tsv', sep='\t', index=False)
    
    print("✓ Created sample_data.tsv and sample_reference.tsv")


def example_level1_simple():
    """Level 1: Simple function interface for 80% of use cases."""
    print("\n" + "="*60)
    print("LEVEL 1: Simple Function Interface")
    print("="*60)
    
    # Simplest possible usage
    params = optimize_filters(
        "sample_data.tsv",
        "sample_reference.tsv",
        n_trials=50,  # Quick for demo
        preset="fast",
        seed=42
    )
    
    print("\nOptimized parameters:")
    for var_type, var_params in params.items():
        print(f"\n{var_type}:")
        for param, value in var_params.items():
            if isinstance(value, float):
                print(f"  {param}: {value:.4f}")
            else:
                print(f"  {param}: {value}")


def example_level2_config():
    """Level 2: Configuration-based for advanced users."""
    print("\n" + "="*60)
    print("LEVEL 2: Configuration Object")
    print("="*60)
    
    # Custom configuration for noisy data
    config = PipelineConfig(
        stage1={'n_trials': 30},  # Quick warmup
        stage2={'min_dnm_count': 5, 'max_dnm_count': 300},  # Outlier removal
        stage3={'n_trials': 50, 'sampler': 'tpe', 'pruner': 'successive_halving'},
        max_workers=4,
        seed=42
    )
    
    print("Configuration:")
    print(f"  Stage 1: {config.stage1.n_trials} warmup trials")
    print(f"  Stage 2: Remove < {config.stage2.min_dnm_count} or > {config.stage2.max_dnm_count} DNMs")
    print(f"  Stage 3: {config.stage3.n_trials} optimization trials with {config.stage3.sampler}")
    
    # Run with configuration
    result = run_optimization(
        "sample_data.tsv",
        "sample_reference.tsv",
        config=config
    )
    
    print(f"\n{result.summary}")


def example_level3_direct():
    """Level 3: Direct stage control for researchers."""
    print("\n" + "="*60)
    print("LEVEL 3: Direct Stage Control")
    print("="*60)
    
    # Load data using abstraction
    data = VariantDataset.from_tsv(Path("sample_data.tsv"))
    reference = VariantDataset.from_tsv(Path("sample_reference.tsv"))
    
    print(f"Loaded {len(data)} variants")
    print(f"Data summary: {data.summary}")
    
    # Calculate targets
    from variant_optimizer.api import calculate_regression_targets
    targets = calculate_regression_targets(reference)
    print(f"\nCalculated targets for {len(targets)} variant types")
    
    # Direct control over stages
    runner = StageRunner()
    
    # Stage 1: Custom warmup
    print("\nStage 1: Warmup...")
    warmup_params = runner.warmup(data, targets, n_trials=20)
    
    # Stage 2: Custom outlier removal
    print("\nStage 2: Outlier removal...")
    clean_data = runner.remove_outliers(
        data, warmup_params,
        min_dnm=10, max_dnm=200
    )
    
    # Stage 3: Custom optimization
    print("\nStage 3: Full optimization...")
    final_params = runner.optimize(clean_data, targets, n_trials=30)
    
    print("\nFinal parameters:")
    for var_type, params in final_params.items():
        print(f"  {var_type}: {len(params)} parameters optimized")


def example_performance_optimizations():
    """Demonstrate performance optimizations."""
    print("\n" + "="*60)
    print("PERFORMANCE OPTIMIZATIONS")
    print("="*60)
    
    # Load data with automatic dtype optimization
    data = VariantDataset.from_tsv(Path("sample_data.tsv"))
    
    # Show memory optimization
    print("\nMemory optimization through categorical dtypes:")
    for col in ['CHROM', 'REF', 'ALT']:
        if col in data.variants.columns:
            dtype = data.variants[col].dtype
            print(f"  {col}: {dtype}")
    
    # Configuration with caching enabled
    config = PipelineConfig(
        use_cache=True,
        cache_dir=Path("./cache"),
        max_workers=4,
        stage3={'pruner': 'successive_halving'}  # Pruning for speed
    )
    
    print(f"\nPerformance settings:")
    print(f"  Caching: {config.use_cache}")
    print(f"  Parallel workers: {config.max_workers}")
    print(f"  Pruning: {config.stage3.pruner}")


def example_config_management():
    """Demonstrate configuration management."""
    print("\n" + "="*60)
    print("CONFIGURATION MANAGEMENT")
    print("="*60)
    
    # Save configuration to file
    config = PipelineConfig.from_preset("noisy_data")
    config.to_yaml(Path("example_config.yaml"))
    print("✓ Saved configuration to example_config.yaml")
    
    # Load and modify
    loaded = PipelineConfig.from_yaml(Path("example_config.yaml"))
    loaded.stage3.n_trials = 1000
    print(f"✓ Loaded and modified config: {loaded.stage3.n_trials} trials")
    
    # Show available presets
    from variant_optimizer.config import PRESETS
    print("\nAvailable presets:")
    for name in PRESETS:
        print(f"  - {name}")


def main():
    """Run all examples."""
    print("VARIANT OPTIMIZER - Refactored Implementation Examples")
    print("="*60)
    
    # Create sample data
    create_sample_data()
    
    # Run examples
    example_level1_simple()
    example_level2_config()
    example_level3_direct()
    example_performance_optimizations()
    example_config_management()
    
    print("\n" + "="*60)
    print("EXAMPLES COMPLETE")
    print("="*60)
    print("\nKey improvements demonstrated:")
    print("✓ Progressive API (3 levels of complexity)")
    print("✓ Performance optimizations (vectorization, caching, pruning)")
    print("✓ Configuration management (YAML, presets, validation)")
    print("✓ Data abstractions (VariantDataset)")
    print("✓ Reproducibility (deterministic seeds)")
    
    # Clean up
    import os
    for f in ['sample_data.tsv', 'sample_reference.tsv', 'example_config.yaml']:
        if os.path.exists(f):
            os.remove(f)


if __name__ == "__main__":
    main()
