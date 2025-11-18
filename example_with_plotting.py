#!/usr/bin/env python3
"""
Example: Running DNM optimization with automatic plotting.

This example demonstrates how to use the dnm-harmoniser package with
automatic plot generation for paternal and maternal age separately.
"""

from dnm_harmoniser import run_optimisation, PipelineConfig
from pathlib import Path

def main():
    """Run optimization with automatic plotting."""

    # Configuration file path
    config_file = Path("/Users/kartikchundru/dnms/ukb/filter.yaml")

    # Data file paths (update these to your actual data paths)
    data_path = "/Users/kartikchundru/dnms/ukb/ukbb_data.tsv"
    reference_path = "/Users/kartikchundru/dnms/decode/decode_reference.tsv"

    # Output directory for results
    output_dir = Path("/Users/kartikchundru/dnms/results")

    print("="*80)
    print("DNM HARMONISER - Optimization with Automatic Plotting")
    print("="*80)
    print(f"\nConfiguration: {config_file}")
    print(f"Data: {data_path}")
    print(f"Reference: {reference_path}")
    print(f"Output: {output_dir}")
    print()

    # Run optimization with automatic plotting
    print("Starting optimization...")
    result = run_optimisation(
        data=data_path,
        reference=reference_path,
        config_file=config_file,
        output_dir=output_dir,
        generate_plots=True  # This is the default
    )

    print("\n" + "="*80)
    print("OPTIMIZATION COMPLETE")
    print("="*80)
    print(result.summary)

    print("\n" + "="*80)
    print("OUTPUT FILES")
    print("="*80)
    print(f"\nResults saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  ✓ optimization_results.png    - 2x3 grid plot (paternal/maternal × SNV/Ins/Del)")
    print("  ✓ filtered_variants.tsv       - Filtered variants passing optimal thresholds")
    print("  ✓ filter_summary.txt          - Summary statistics and retention rates")
    print("  ✓ best_parameters.txt         - Human-readable optimal parameters")
    print("  ✓ optimal_params.yaml         - Machine-readable optimal parameters")
    print("  ✓ summary.txt                 - Complete optimization summary")

    print("\n" + "="*80)
    print("BEST PARAMETERS BY VARIANT TYPE")
    print("="*80)
    for var_type, params in result.best_params.items():
        print(f"\n{var_type}:")
        print("-" * 40)
        for param, value in sorted(params.items()):
            if isinstance(value, float):
                print(f"  {param:25s} = {value:.4f}")
            else:
                print(f"  {param:25s} = {value}")

    print("\n✓ Done! Check the output directory for plots and filtered variants.\n")

    return result


if __name__ == "__main__":
    # Simple usage example
    print("""
SIMPLE USAGE:
=============

from dnm_harmoniser import optimize_filters

params = optimize_filters(
    data_path="ukbb_data.tsv",
    reference_path="decode_reference.tsv",
    output_dir="results/",
    preset="balanced",
    n_trials=100
)

This will automatically:
  1. Run optimization
  2. Generate 2x3 grid plots (paternal/maternal × SNV/Insertion/Deletion)
  3. Save filtered variants to TSV
  4. Save summary statistics
  5. Save best parameters

Press Ctrl+C to skip, or wait to run the full example...
    """)

    import time
    try:
        time.sleep(3)
        main()
    except KeyboardInterrupt:
        print("\n\nExample skipped.")
    except FileNotFoundError as e:
        print(f"\n\nError: Could not find file: {e}")
        print("\nPlease update the file paths in this script to match your data locations.")
    except Exception as e:
        print(f"\n\nError occurred: {e}")
        print("\nMake sure you have:")
        print("  1. Updated file paths to your actual data")
        print("  2. Installed the package: pip install -e .")
        print("  3. Your data files have the required columns")
