#!/usr/bin/env python3
"""
Example: Three-Stage Optimization for Noisy Data
================================================

This script demonstrates how to use the warmup-based outlier removal
to handle very noisy genomic data.

The three stages:
1. Fast warmup optimization (50 trials)
2. Outlier removal based on FILTERED DNM counts
3. Full optimization on cleaned data (500 trials)
"""

from optimise_optuna_simplified import main

# ============================================================================
# EXAMPLE 1: Your Exact Use Case
# ============================================================================
# Remove individuals with <5 or >300 FILTERED DNMs across all variant types
# Uses warmup to determine filters BEFORE counting DNMs

print("="*70)
print("EXAMPLE 1: Three-Stage Optimization for Noisy Data")
print("="*70)
print()
print("Settings:")
print("  - Warmup: 50 trials (fast)")
print("  - Full optimization: 500 trials")
print("  - Remove if < 5 FILTERED DNMs")
print("  - Remove if > 300 FILTERED DNMs")
print()

results = main(
    data_path='data_snv_dnms_for_filtering.tsv',
    reference_path='reference_parages.txt',
    
    # Outlier removal thresholds (based on FILTERED DNM counts)
    min_dnm_count=5,        # Remove individuals with < 5 filtered DNMs
    max_dnm_count=300,      # Remove individuals with > 300 filtered DNMs
    
    # Optimization settings
    warmup_trials=50,       # Fast warmup to get initial filters
    n_trials=500,          # Full optimization after outlier removal
    
    random_seed=2025
)

print("\n" + "="*70)
print("RESULTS")
print("="*70)
for var_type, params in results.items():
    print(f"\n{var_type} optimal parameters:")
    for param_name, param_value in params.items():
        if isinstance(param_value, float):
            print(f"  {param_name}: {param_value:.4f}")
        else:
            print(f"  {param_name}: {param_value}")


# ============================================================================
# EXAMPLE 2: Conservative Outlier Removal
# ============================================================================
# Only remove extreme outliers, keep more individuals

print("\n\n" + "="*70)
print("EXAMPLE 2: Conservative Outlier Removal")
print("="*70)
print()
print("Settings:")
print("  - Less aggressive outlier removal")
print("  - Remove if < 3 FILTERED DNMs (very low)")  
print("  - Remove if > 500 FILTERED DNMs (very high)")
print()

results_conservative = main(
    data_path='data_snv_dnms_for_filtering.tsv',
    reference_path='reference_parages.txt',
    min_dnm_count=3,        # More lenient lower bound
    max_dnm_count=500,      # More lenient upper bound
    warmup_trials=50,
    n_trials=500,
    random_seed=2025
)


# ============================================================================
# EXAMPLE 3: Aggressive Cleaning
# ============================================================================
# Remove more individuals for a very clean dataset

print("\n\n" + "="*70)
print("EXAMPLE 3: Aggressive Outlier Removal")
print("="*70)
print()
print("Settings:")
print("  - Aggressive outlier removal for cleanest data")
print("  - Remove if < 10 FILTERED DNMs")
print("  - Remove if < 200 FILTERED DNMs")
print()

results_aggressive = main(
    data_path='data_snv_dnms_for_filtering.tsv',
    reference_path='reference_parages.txt',
    min_dnm_count=10,       # Higher lower bound
    max_dnm_count=200,      # Lower upper bound
    warmup_trials=50,
    n_trials=500,
    random_seed=2025
)


# ============================================================================
# EXAMPLE 4: Fast Warmup, More Thorough Full Optimization
# ============================================================================
# Quick warmup but extensive full optimization

print("\n\n" + "="*70)
print("EXAMPLE 4: Quick Warmup, Thorough Final Optimization")
print("="*70)
print()
print("Settings:")
print("  - Very fast warmup: 30 trials")
print("  - Very thorough full optimization: 1000 trials")
print()

results_thorough = main(
    data_path='data_snv_dnms_for_filtering.tsv',
    reference_path='reference_parages.txt',
    min_dnm_count=5,
    max_dnm_count=300,
    warmup_trials=30,       # Even faster warmup
    n_trials=1000,          # More thorough full optimization
    random_seed=2025
)


# ============================================================================
# EXAMPLE 5: Custom Sample Column Names
# ============================================================================
# When your data uses different column names for sample IDs

print("\n\n" + "="*70)
print("EXAMPLE 5: Custom Sample Column Names")
print("="*70)
print()
print("Settings:")
print("  - Data file uses 'patient_id' for samples")
print("  - Reference file uses 'subject_id' for samples")
print()

results_custom_cols = main(
    data_path='data_snv_dnms_for_filtering.tsv',
    reference_path='reference_parages.txt',
    
    # Custom sample column names
    sample_col_data='patient_id',       # Column name in data file
    sample_col_reference='subject_id',  # Column name in reference file
    
    min_dnm_count=5,
    max_dnm_count=300,
    warmup_trials=50,
    n_trials=500,
    random_seed=2025
)


# ============================================================================
# EXAMPLE 6: No Warmup - Just Remove Outliers Based on Raw Counts
# ============================================================================
# Sometimes you might want to skip warmup and just remove based on raw counts
# (Not recommended for noisy data, but included for completeness)

print("\n\n" + "="*70)
print("EXAMPLE 6: No Warmup (Not Recommended for Noisy Data)")
print("="*70)
print()
print("Settings:")
print("  - No warmup optimization")
print("  - Outlier removal based on RAW DNM counts")
print("  - WARNING: Less robust for noisy data!")
print()

results_no_warmup = main(
    data_path='data_snv_dnms_for_filtering.tsv',
    reference_path='reference_parages.txt',
    min_dnm_count=5,
    max_dnm_count=300,
    warmup_trials=None,     # No warmup! Will use raw counts
    n_trials=500,
    random_seed=2025
)


# ============================================================================
# EXAMPLE 7: Only Warmup, No Outlier Removal
# ============================================================================
# Use warmup to speed up optimization, but don't remove outliers

print("\n\n" + "="*70)
print("EXAMPLE 7: Warmup for Speed, But No Outlier Removal")
print("="*70)
print()
print("Settings:")
print("  - Fast warmup: 50 trials")
print("  - Full optimization: 500 trials")
print("  - No outlier removal")
print("  - Good for high-quality data where you just want faster optimization")
print()

results_warmup_only = main(
    data_path='data_snv_dnms_for_filtering.tsv',
    reference_path='reference_parages.txt',
    warmup_trials=50,       # Warmup enabled
    n_trials=500,
    min_dnm_count=None,     # No outlier removal
    max_dnm_count=None,     # No outlier removal
    random_seed=2025
)
