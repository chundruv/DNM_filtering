#!/usr/bin/env python3
"""Debug optimization issues."""

import pandas as pd
from dnm_harmoniser import PipelineConfig, VariantDataset
import optuna

# Load configuration
config = PipelineConfig.from_yaml('/Users/kartikchundru/dnms/ukb/filter.yaml')

print("="*80)
print("CONFIGURATION CHECK")
print("="*80)

# Check optimization columns
opt_cols = config.optimisation.get_optimisation_columns()
print(f"\nOptimization columns configured: {len(opt_cols)}")
for col in opt_cols:
    constraint_info = ""
    if col.range_constraint:
        if hasattr(col.range_constraint, 'lower'):  # SymmetricRangeConstraint
            constraint_info = f" (symmetric: {col.range_constraint.min}-{col.range_constraint.max}, scale={col.range_constraint.scale})"
        else:  # RangeConstraint
            constraint_info = f" (range: {col.range_constraint.min}-{col.range_constraint.max})"

    linked_info = f" [linked to: {col.linked_to}]" if col.linked_to else ""
    variant_info = f" [variants: {col.variant_types}]" if col.variant_types else ""

    print(f"  {col.name}: {col.optimisation} ({col.dtype}){constraint_info}{linked_info}{variant_info}")

print("\n" + "="*80)
print("TESTING PARAMETER SUGGESTION")
print("="*80)

# Create a mock trial to test parameter suggestion
class MockTrial:
    def __init__(self):
        self.params = {}

    def suggest_int(self, name, low, high):
        mid = (low + high) // 2
        self.params[name] = mid
        print(f"  suggest_int('{name}', {low}, {high}) → {mid}")
        return mid

    def suggest_float(self, name, low, high):
        mid = (low + high) / 2
        self.params[name] = mid
        print(f"  suggest_float('{name}', {low:.2f}, {high:.2f}) → {mid:.2f}")
        return mid

# Test parameter suggestion (requires data - let's create a simple mock)
print("\nNote: To test actual parameter suggestion, we need to load data.")
print("Parameters will be in format:")
print("  - For 'minimum' optimization: max_{column_name}")
print("  - For 'maximum' optimization: min_{column_name}")
print("  - For 'range' optimization: min_{column_name} and max_{column_name}")

print("\n" + "="*80)
print("EXPECTED PARAMETER NAMES")
print("="*80)

for col in opt_cols:
    if col.optimisation == 'minimum':
        print(f"  {col.name} → max_{col.name}")
    elif col.optimisation == 'maximum':
        print(f"  {col.name} → min_{col.name}")
    elif col.optimisation == 'range':
        print(f"  {col.name} → min_{col.name} AND max_{col.name}")

print("\n" + "="*80)
print("RECOMMENDATIONS")
print("="*80)
print("""
1. Update your plotting code to use new parameter names:

   OLD: best_params['min_cnn_prob']
   NEW: best_params['max_DeNovoCNN_prob']

   OLD: best_params['min_dnm']
   NEW: best_params['max_DNM']

   OLD: best_params['min_mq']
   NEW: best_params['max_MQ']

   OLD: best_params['max_nparaadn0']
   NEW: best_params['max_nparAADn0']

   OLD: best_params['min_cov'] (for all coverage columns)
   NEW: best_params['min_child_coverage']
        best_params['min_father_coverage']
        best_params['min_mother_coverage']

   OLD: best_params['min_vaf']
   NEW: best_params['min_VAF'] and best_params['max_VAF']

   OLD: best_params['min_imf']
   NEW: best_params['min_IMF'] and best_params['max_IMF']

2. Note: Linked columns (father_coverage, mother_coverage) will share values

3. Symmetric range constraints:
   - VAF: lower=25, scale=100 → min_VAF and max_VAF will be suggested
   - You don't need to calculate 100-min_vaf anymore - use max_VAF directly
""")
