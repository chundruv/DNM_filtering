#!/usr/bin/env python3
"""Test that all three stages are actually called."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from dnm_harmoniser import PipelineConfig, run_optimisation
import logging

# Set up logging to see everything
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s - %(name)s - %(message)s'
)

# Load your config
config_path = "/Users/kartikchundru/dnms/ukb/filter.yaml"
config = PipelineConfig.from_yaml(config_path)

print("="*80)
print("CONFIG CHECK")
print("="*80)
print(f"Stage 1 enabled: {config.stage1.enabled}, trials: {config.stage1.n_trials}")
print(f"Stage 2 enabled: {config.stage2.enabled}, min/max DNM: {config.stage2.min_dnm_count}/{config.stage2.max_dnm_count}")
print(f"Stage 3 enabled: {config.stage3.enabled}, trials: {config.stage3.n_trials}")

print("\n" + "="*80)
print("RUNNING OPTIMIZATION (with 10 trials for testing)")
print("="*80)

# Override to use fewer trials for testing
config.stage1.n_trials = 10
config.stage3.n_trials = 10

try:
    result = run_optimisation(
        data="/Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv",
        reference="/Users/kartikchundru/dnms/decode_parages.txt",
        config=config,
        output_dir="test_output",
        generate_plots=False  # Disable plots for quick test
    )

    print("\n" + "="*80)
    print("RESULT")
    print("="*80)
    print(f"Best parameters found for {len(result.best_params)} variant types")
    print(f"Warmup parameters: {result.warmup_params is not None}")
    print(f"Outliers removed: {result.n_individuals_removed}")

except Exception as e:
    print(f"\nError: {e}")
    import traceback
    traceback.print_exc()
