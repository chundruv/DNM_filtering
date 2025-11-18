#!/usr/bin/env python3
"""Diagnostic script to debug infinite trial values."""

import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from dnm_harmoniser import PipelineConfig, VariantDataset
import numpy as np

def diagnose_issue(data_path, reference_path, config_path):
    """Diagnose why optimization returns inf."""

    print("="*80)
    print("DIAGNOSTIC REPORT: Debugging Infinite Trial Values")
    print("="*80)

    # Load configuration
    print("\n1. LOADING CONFIGURATION")
    print("-"*80)
    try:
        config = PipelineConfig.from_yaml(config_path)
        print(f"✓ Config loaded from {config_path}")

        opt_cols = config.optimisation.get_optimisation_columns()
        print(f"✓ Found {len(opt_cols)} optimization columns:")
        for col in opt_cols:
            print(f"    - {col.name}: {col.optimisation} ({col.dtype})")

    except Exception as e:
        print(f"✗ Failed to load config: {e}")
        return

    # Load data
    print("\n2. LOADING DATA")
    print("-"*80)
    try:
        print(f"Loading {data_path}...")
        data = VariantDataset.from_tsv(
            Path(data_path),
            sample_col=config.optimisation.sample_id_column_data,
            paternal_age_col=config.optimisation.paternal_age_column_data,
            maternal_age_col=config.optimisation.maternal_age_column_data,
            reference_col=config.optimisation.reference_column_data,
            alternate_col=config.optimisation.alternate_column_data
        )
        print(f"✓ Loaded {len(data)} variants")
        print(f"  Columns: {list(data.variants.columns)}")
        print(f"  Variant types: {data.summary['variant_types']}")

    except Exception as e:
        print(f"✗ Failed to load data: {e}")
        import traceback
        traceback.print_exc()
        return

    # Load reference
    print(f"\nLoading {reference_path}...")
    try:
        reference = VariantDataset.from_tsv(
            Path(reference_path),
            sample_col=config.optimisation.sample_id_column_reference,
            paternal_age_col=config.optimisation.paternal_age_column_reference,
            maternal_age_col=config.optimisation.maternal_age_column_reference,
            reference_col=config.optimisation.reference_column_reference,
            alternate_col=config.optimisation.alternate_column_reference
        )
        print(f"✓ Loaded {len(reference)} reference variants")
        print(f"  Variant types: {reference.summary['variant_types']}")

    except Exception as e:
        print(f"✗ Failed to load reference: {e}")
        import traceback
        traceback.print_exc()
        return

    # Check columns exist in data
    print("\n3. CHECKING COLUMN AVAILABILITY")
    print("-"*80)
    missing_cols = []
    for col_config in opt_cols:
        col = col_config.name
        if col in data.variants.columns:
            dtype = data.variants[col].dtype
            non_null = data.variants[col].notna().sum()
            print(f"✓ {col:25s} - dtype: {dtype}, non-null: {non_null}/{len(data)}")
        else:
            print(f"✗ {col:25s} - MISSING!")
            missing_cols.append(col)

    if missing_cols:
        print(f"\n⚠ WARNING: Missing columns: {missing_cols}")
        print("  These columns are in your config but not in your data!")
        print("  Check column names in your data files.")

    # Check required metadata columns
    print("\n4. CHECKING REQUIRED METADATA COLUMNS")
    print("-"*80)
    required = ['SAMPLE', 'paternal_age', 'maternal_age', 'REF', 'ALT', 'var_type']
    for col in required:
        if col in data.variants.columns:
            if col == 'var_type':
                types = data.variants[col].unique()
                print(f"✓ {col:25s} - values: {types}")
            elif col in ['SAMPLE', 'REF', 'ALT']:
                unique = data.variants[col].nunique()
                print(f"✓ {col:25s} - {unique} unique values")
            else:
                min_val = data.variants[col].min()
                max_val = data.variants[col].max()
                print(f"✓ {col:25s} - range: {min_val:.1f} to {max_val:.1f}")
        else:
            print(f"✗ {col:25s} - MISSING!")

    # Check data ranges for optimization columns
    print("\n5. DATA RANGES FOR OPTIMIZATION COLUMNS")
    print("-"*80)
    for col_config in opt_cols:
        col = col_config.name
        if col not in data.variants.columns:
            continue

        col_data = pd.to_numeric(data.variants[col], errors='coerce').dropna()
        if len(col_data) == 0:
            print(f"✗ {col:25s} - All values are NULL or non-numeric!")
            continue

        print(f"\n{col} ({col_config.optimisation} optimization):")
        print(f"  Count:  {len(col_data)}")
        print(f"  Min:    {col_data.min()}")
        print(f"  25%:    {col_data.quantile(0.25)}")
        print(f"  Median: {col_data.median()}")
        print(f"  75%:    {col_data.quantile(0.75)}")
        print(f"  Max:    {col_data.max()}")
        print(f"  NaN:    {data.variants[col].isna().sum()}")

    # Test filtering with dummy parameters
    print("\n6. TESTING SAMPLE FILTERING")
    print("-"*80)

    for var_type in ['SNV', 'Insertion', 'Deletion']:
        var_data = data.filter_by_type(var_type)
        if len(var_data) == 0:
            print(f"\n{var_type}: No variants of this type")
            continue

        print(f"\n{var_type}: {len(var_data)} variants")

        # Create sample parameters (mid-range values)
        test_params = {}
        for col_config in opt_cols:
            col = col_config.name
            if col not in var_data.variants.columns:
                continue

            # Skip variant-specific columns if they don't apply
            if col_config.variant_types and var_type not in col_config.variant_types:
                continue

            col_data = pd.to_numeric(var_data.variants[col], errors='coerce').dropna()
            if len(col_data) == 0:
                continue

            if col_config.optimisation == 'minimum':
                # Use median as test threshold
                test_params[f'max_{col}'] = float(col_data.median())
            elif col_config.optimisation == 'maximum':
                # Use median as test threshold
                test_params[f'min_{col}'] = float(col_data.median())
            elif col_config.optimisation == 'range':
                # Use 25th and 75th percentiles
                test_params[f'min_{col}'] = float(col_data.quantile(0.25))
                test_params[f'max_{col}'] = float(col_data.quantile(0.75))

        print(f"  Test parameters:")
        for param, value in sorted(test_params.items()):
            print(f"    {param}: {value:.4f}")

        # Apply filters
        filtered = var_data.apply_filters(test_params)
        retention = len(filtered) / len(var_data) * 100 if len(var_data) > 0 else 0
        print(f"  Filtered: {len(var_data)} → {len(filtered)} ({retention:.1f}% retained)")

        if len(filtered) == 0:
            print(f"  ✗ WARNING: All variants filtered out!")
            print(f"     This will cause inf values in optimization")
            print(f"     Check your parameter ranges")
        else:
            # Try to count per sample
            if 'SAMPLE' in filtered.variants.columns:
                sample_counts = filtered.count_by_sample()
                print(f"  Sample count range: {sample_counts.min()}-{sample_counts.max()}")
                print(f"  Mean DNMs per sample: {sample_counts.mean():.1f}")

    # Check if reference data has required columns for regression
    print("\n7. CHECKING REFERENCE DATA FOR REGRESSION")
    print("-"*80)
    ref_required = ['SAMPLE', 'paternal_age', 'maternal_age', 'var_type']
    ref_ok = True
    for col in ref_required:
        if col in reference.variants.columns:
            print(f"✓ {col:25s} - present")
        else:
            print(f"✗ {col:25s} - MISSING!")
            ref_ok = False

    if ref_ok:
        print("\n  Reference data by variant type:")
        for var_type in ['SNV', 'Insertion', 'Deletion']:
            ref_subset = reference.filter_by_type(var_type)
            if len(ref_subset) > 0:
                n_samples = ref_subset.variants['SAMPLE'].nunique()
                age_range = (ref_subset.variants['paternal_age'].min(),
                           ref_subset.variants['paternal_age'].max())
                print(f"    {var_type}: {len(ref_subset)} variants, {n_samples} samples, "
                      f"age range {age_range[0]:.0f}-{age_range[1]:.0f}")

    # Summary and recommendations
    print("\n" + "="*80)
    print("SUMMARY & RECOMMENDATIONS")
    print("="*80)

    issues = []

    if missing_cols:
        issues.append(f"Missing columns in data: {missing_cols}")
        print(f"\n✗ CRITICAL: Your data is missing optimization columns: {missing_cols}")
        print(f"  → Check that column names in filter.yaml match your data file")

    if not ref_ok:
        issues.append("Missing columns in reference data")
        print(f"\n✗ CRITICAL: Reference data is missing required columns")
        print(f"  → Reference must have: SAMPLE, paternal_age, maternal_age, var_type")

    # Check if any test filtering resulted in 0 variants
    print("\nTo fix the 'inf' issue:")
    print("  1. Verify column names match between config and data files")
    print("  2. Ensure all optimization columns have numeric values")
    print("  3. Check that filter ranges allow some variants to pass")
    print("  4. Verify data has 'var_type' column with values: SNV, Insertion, Deletion")
    print("  5. Run with verbose flag: dnm-harmoniser run ... -vv")

    if not issues:
        print("\n✓ No critical issues found!")
        print("  If you're still getting inf, try running with -vv for detailed logs")

    print()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python diagnose_inf_issue.py <data.tsv> <reference.tsv> <config.yaml>")
        print("\nExample:")
        print("  python diagnose_inf_issue.py \\")
        print("    /Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv \\")
        print("    /Users/kartikchundru/dnms/decode_parages.txt \\")
        print("    /Users/kartikchundru/dnms/ukb/filter.yaml")
        sys.exit(1)

    data_path = sys.argv[1]
    reference_path = sys.argv[2]
    config_path = sys.argv[3]

    diagnose_issue(data_path, reference_path, config_path)
