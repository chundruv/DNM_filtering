#!/usr/bin/env python3
"""Deep debugging of optimization pipeline."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

import pandas as pd
import numpy as np
from dnm_harmoniser import PipelineConfig, VariantDataset

def deep_debug():
    """Trace through the entire optimization process."""

    print("="*80)
    print("DEEP DEBUG: Tracing Optimization Pipeline")
    print("="*80)

    # Paths
    data_path = "/Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv"
    reference_path = "/Users/kartikchundru/dnms/decode_parages.txt"
    config_path = "/Users/kartikchundru/dnms/ukb/filter.yaml"

    # Load config
    print("\n1. LOADING CONFIGURATION")
    print("-"*80)
    config = PipelineConfig.from_yaml(config_path)
    opt_cols = config.optimisation.get_optimisation_columns()
    print(f"Optimization columns: {[c.name for c in opt_cols]}")

    # Check linked columns
    linked_groups = config.optimisation.get_linked_column_groups()
    print(f"Linked column groups: {[[c.name for c in group] for group in linked_groups]}")

    # Load data
    print("\n2. LOADING DATA")
    print("-"*80)
    data = VariantDataset.from_tsv(
        Path(data_path),
        sample_col=config.optimisation.sample_id_column_data,
        paternal_age_col=config.optimisation.paternal_age_column_data,
        maternal_age_col=config.optimisation.maternal_age_column_data,
        reference_col=config.optimisation.reference_column_data,
        alternate_col=config.optimisation.alternate_column_data
    )
    print(f"Loaded {len(data)} variants")

    # Test parameter suggestion for SNVs
    print("\n3. TESTING PARAMETER SUGGESTION FOR SNVs")
    print("-"*80)
    snv_data = data.filter_by_type('SNV')
    print(f"SNV data: {len(snv_data)} variants")

    # Simulate what _suggest_params does
    print("\nSimulating parameter suggestion:")
    test_params = {}

    for col_config in opt_cols:
        col = col_config.name
        if col not in snv_data.variants.columns:
            print(f"  Skipping {col} - not in data")
            continue

        # Skip variant-specific columns
        if col_config.variant_types and 'SNV' not in col_config.variant_types:
            print(f"  Skipping {col} - not for SNVs")
            continue

        col_data = pd.to_numeric(snv_data.variants[col], errors='coerce').dropna()
        if len(col_data) == 0:
            print(f"  Skipping {col} - no valid data")
            continue

        print(f"\n  {col} ({col_config.optimisation}):")
        print(f"    Data type: {col_data.dtype}")
        print(f"    Min: {col_data.min()}")
        print(f"    10%: {col_data.quantile(0.1)}")
        print(f"    Median: {col_data.median()}")
        print(f"    90%: {col_data.quantile(0.9)}")
        print(f"    Max: {col_data.max()}")

        if col_config.optimisation == 'minimum':
            lower = col_data.min()
            upper = col_data.quantile(0.9)
            print(f"    → Will suggest max_{col} in range [{lower}, {upper}]")
            test_params[f'max_{col}'] = col_data.quantile(0.9)

        elif col_config.optimisation == 'maximum':
            lower = col_data.quantile(0.1)
            upper = col_data.max()
            print(f"    → Will suggest min_{col} in range [{lower}, {upper}]")
            test_params[f'min_{col}'] = col_data.quantile(0.1)

        elif col_config.optimisation == 'range':
            if col_config.range_constraint:
                lower_bound = col_config.range_constraint.min
                upper_bound = col_config.range_constraint.max
                print(f"    → Range constraint: [{lower_bound}, {upper_bound}]")
            else:
                lower_bound = col_data.min()
                upper_bound = col_data.max()
                print(f"    → Will suggest range in [{lower_bound}, {upper_bound}]")
            test_params[f'min_{col}'] = lower_bound if col_config.range_constraint else col_data.quantile(0.25)
            test_params[f'max_{col}'] = upper_bound if col_config.range_constraint else col_data.quantile(0.75)

    # Apply linked columns logic
    print("\n4. APPLYING LINKED COLUMNS LOGIC")
    print("-"*80)
    for group in linked_groups:
        if len(group) < 2:
            continue
        print(f"\nLinked group: {[c.name for c in group]}")
        first_col = group[0].name
        param_keys = [k for k in test_params.keys() if first_col in k]
        print(f"  First column params: {param_keys}")

        if param_keys:
            for key in param_keys:
                param_value = test_params[key]
                prefix = key.split('_')[0]
                print(f"  {key} = {param_value}")

                for linked_col in group[1:]:
                    linked_key = f'{prefix}_{linked_col.name}'
                    if linked_key in test_params:
                        old_value = test_params[linked_key]
                        test_params[linked_key] = param_value
                        print(f"    {linked_key}: {old_value} → {param_value}")

    print("\n5. FINAL TEST PARAMETERS")
    print("-"*80)
    for param, value in sorted(test_params.items()):
        print(f"  {param}: {value}")

    # Test filtering
    print("\n6. TESTING FILTERING")
    print("-"*80)
    print(f"Before filtering: {len(snv_data)} variants")

    # Apply filters step by step
    filtered = snv_data.variants.copy()
    for param, value in sorted(test_params.items()):
        before = len(filtered)

        if param.startswith('min_'):
            col = param[4:]
            if col in filtered.columns:
                filtered = filtered[filtered[col] >= value]
                print(f"  After {param} >= {value}: {len(filtered)} variants ({before - len(filtered)} removed)")

        elif param.startswith('max_'):
            col = param[4:]
            if col in filtered.columns:
                filtered = filtered[filtered[col] <= value]
                print(f"  After {param} <= {value}: {len(filtered)} variants ({before - len(filtered)} removed)")

    print(f"\nFinal: {len(filtered)} variants ({len(filtered)/len(snv_data)*100:.1f}% retained)")

    if len(filtered) == 0:
        print("\n⚠ ALL VARIANTS FILTERED OUT - This will cause inf!")
        print("\nMost restrictive filters:")

        # Re-run to find most restrictive
        filtered_test = snv_data.variants.copy()
        removals = []
        for param, value in sorted(test_params.items()):
            before = len(filtered_test)
            if param.startswith('min_'):
                col = param[4:]
                if col in filtered_test.columns:
                    filtered_test = filtered_test[filtered_test[col] >= value]
            elif param.startswith('max_'):
                col = param[4:]
                if col in filtered_test.columns:
                    filtered_test = filtered_test[filtered_test[col] <= value]
            removed = before - len(filtered_test)
            if removed > 0:
                removals.append((param, value, removed, removed/before*100))

        removals.sort(key=lambda x: x[2], reverse=True)
        for param, value, removed, pct in removals[:5]:
            print(f"  {param} = {value}: removed {removed} variants ({pct:.1f}%)")

    # Test regression
    if len(filtered) > 0:
        print("\n7. TESTING REGRESSION")
        print("-"*80)
        filtered_df = pd.DataFrame(filtered)

        # Count DNMs per sample
        if 'SAMPLE' in filtered_df.columns:
            dnm_counts = filtered_df.groupby('SAMPLE').size().rename('dnm_count')
            print(f"DNM counts: min={dnm_counts.min()}, max={dnm_counts.max()}, mean={dnm_counts.mean():.1f}")

            # Get parental ages
            if 'paternal_age' in filtered_df.columns and 'maternal_age' in filtered_df.columns:
                ages = filtered_df[['SAMPLE', 'paternal_age', 'maternal_age']].drop_duplicates().set_index('SAMPLE')
                regression_data = ages.join(dnm_counts, how='left').fillna(0)

                print(f"Regression data: {len(regression_data)} samples")
                print(f"  DNM range: {regression_data['dnm_count'].min()}-{regression_data['dnm_count'].max()}")
                print(f"  Age range: {regression_data['paternal_age'].min()}-{regression_data['paternal_age'].max()}")

                # Try regression
                try:
                    import statsmodels.formula.api as smf
                    model = smf.ols('dnm_count ~ paternal_age + maternal_age', data=regression_data).fit()
                    print(f"\nRegression coefficients:")
                    print(f"  Intercept: {model.params['Intercept']:.4f}")
                    print(f"  Paternal age: {model.params['paternal_age']:.4f}")
                    print(f"  Maternal age: {model.params['maternal_age']:.4f}")

                    # Compare to targets
                    print(f"\nTarget coefficients (from reference):")
                    print(f"  Intercept: 11.5293")
                    print(f"  Paternal age: 1.4242")
                    print(f"  Maternal age: 0.3596")

                    # Calculate MSE
                    targets = np.array([11.5293, 1.4242, 0.3596])
                    params = model.params.values
                    mse = np.mean((params - targets) ** 2)
                    print(f"\nMSE: {mse:.6f}")

                except Exception as e:
                    print(f"Regression failed: {e}")

if __name__ == "__main__":
    deep_debug()
