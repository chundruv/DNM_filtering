# Parameter Naming Guide: Fix for `inf` Values and Missing Plots

## Problem

After updating to the new column configuration system, optimization trials are returning `inf` values and no plots are generated because:

1. **Parameter names have changed** - The system now uses exact column names from configuration
2. **Optimization logic inverted** - "Minimum" optimization now uses `max_*` thresholds (counter-intuitive but correct)

## Understanding the New Parameter Names

### Old System (Hardcoded)
```python
# Old parameter names were shortened/hardcoded
best_params = {
    'min_cnn_prob': 0.5,
    'min_dnm': 10.0,
    'min_mq': 40,
    'max_nparaadn0': 5,
    'min_cov': 10,  # Single parameter for all coverage columns
    'min_vaf': 30,  # Only lower bound
}
```

### New System (Configuration-Based)
```python
# New parameter names use exact column names from config
best_params = {
    'max_DeNovoCNN_prob': 0.5,      # Note: max for "minimum" optimization!
    'max_DNM': 10.0,                 # Note: max for "minimum" optimization!
    'max_MQ': 40,                    # Note: max for "minimum" optimization!
    'max_nparAADn0': 5,              # Correct: max threshold
    'min_child_coverage': 10,        # Note: min for "maximum" optimization!
    'min_father_coverage': 8,        # Separate parameter for father
    'min_mother_coverage': 8,        # Separate parameter for mother
    'min_VAF': 30,                   # Lower bound
    'max_VAF': 70,                   # Upper bound (NEW!)
}
```

## Why the Naming Seems Backwards

### Minimum Optimization → Uses `max_*` Threshold

For columns where **lower is better** (like error rates, quality scores):
- Configuration: `optimisation: minimum`
- Generated parameter: `max_{column_name}`
- Logic: Keep values **≤ threshold** (filter OUT high/bad values)

```yaml
# Config
- name: DeNovoCNN_prob
  dtype: float
  optimisation: minimum  # Lower probability is better
```

```python
# Generates: max_DeNovoCNN_prob
# Filtering: df[df['DeNovoCNN_prob'] <= max_DeNovoCNN_prob]
# Keeps: Low probabilities (good variants)
# Filters out: High probabilities (bad variants)
```

### Maximum Optimization → Uses `min_*` Threshold

For columns where **higher is better** (like coverage, depth):
- Configuration: `optimisation: maximum`
- Generated parameter: `min_{column_name}`
- Logic: Keep values **≥ threshold** (filter OUT low/bad values)

```yaml
# Config
- name: child_coverage
  dtype: int
  optimisation: maximum  # Higher coverage is better
```

```python
# Generates: min_child_coverage
# Filtering: df[df['child_coverage'] >= min_child_coverage]
# Keeps: High coverage (good variants)
# Filters out: Low coverage (bad variants)
```

## Complete Parameter Mapping

### Your Configuration → Generated Parameters

| Column | Optimization | Generated Parameter(s) | Filter Logic |
|--------|--------------|------------------------|--------------|
| `DeNovoCNN_prob` | minimum | `max_DeNovoCNN_prob` | `≤ threshold` |
| `DNM` | minimum | `max_DNM` | `≤ threshold` |
| `MQ` | minimum | `max_MQ` | `≤ threshold` |
| `nparAADn0` | minimum | `max_nparAADn0` | `≤ threshold` |
| `child_coverage` | maximum | `min_child_coverage` | `≥ threshold` |
| `father_coverage` | maximum | `min_father_coverage` | `≥ threshold` (linked to mother) |
| `mother_coverage` | maximum | `min_mother_coverage` | `≥ threshold` (linked to father) |
| `VAF` | range | `min_VAF`, `max_VAF` | `≥ min AND ≤ max` |
| `IMF` | range | `min_IMF`, `max_IMF` | `≥ min AND ≤ max` (indels only) |

## Fixing Your Plotting Code

### Old Code Problems

```python
# ❌ OLD - Won't work with new system
filtered_ukbb_df = ukbb_df[
    (ukbb_df['DeNovoCNN_prob'] >= best_params['min_cnn_prob']) &  # Wrong param name
    (ukbb_df['DNM'] >= best_params['min_dnm']) &                   # Wrong param name
    (ukbb_df['VAF'] >= best_params['min_vaf']) &                   # Incomplete - missing max
    (ukbb_df['VAF'] <= 100-best_params['min_vaf'])                 # Manual calculation not needed
]
```

### Updated Code

```python
# ✅ NEW - Correct parameter names and logic
filtered_ukbb_df = ukbb_df[
    (ukbb_df['DeNovoCNN_prob'] <= best_params['max_DeNovoCNN_prob']) &  # max threshold!
    (ukbb_df['DNM'] <= best_params['max_DNM']) &                         # max threshold!
    (ukbb_df['MQ'] <= best_params['max_MQ']) &                           # max threshold!
    (ukbb_df['nparAADn0'] <= best_params['max_nparAADn0']) &
    (ukbb_df['child_coverage'] >= best_params['min_child_coverage']) &
    (ukbb_df['father_coverage'] >= best_params['min_father_coverage']) &
    (ukbb_df['mother_coverage'] >= best_params['min_mother_coverage']) &
    (ukbb_df['VAF'] >= best_params['min_VAF']) &                         # Both bounds
    (ukbb_df['VAF'] <= best_params['max_VAF'])                           # provided!
]
```

## Missing Column: FS (FisherStrand)

Your old code references `FS` but it's not in your configuration. Add it if needed:

```yaml
# Add to filter.yaml
columns:
  - name: FS
    dtype: float
    optimisation: minimum  # Lower Fisher Strand bias is better
```

This will generate `max_FS` parameter.

## Updated Plotting Function

See [updated_plot_function.py](updated_plot_function.py) for a complete, corrected version.

Key changes:
1. Uses exact column names from configuration
2. Applies correct comparison operators
3. Handles symmetric range constraints properly
4. Separate coverage parameters for each sample type
5. Debug output to verify parameters

## Debugging Inf Values

If you're still getting `inf` values, check:

### 1. Data Column Names
```python
import pandas as pd

ukbb_df = pd.read_csv('your_data.tsv', sep='\t')
print("Available columns:", ukbb_df.columns.tolist())

# Ensure these match your config:
required_cols = [
    'SAMPLE', 'DeNovoCNN_prob', 'DNM', 'MQ', 'nparAADn0',
    'child_coverage', 'father_coverage', 'mother_coverage',
    'VAF', 'REF', 'ALT', 'paternal_age', 'maternal_age'
]
missing = [c for c in required_cols if c not in ukbb_df.columns]
if missing:
    print(f"Missing columns: {missing}")
```

### 2. Data Value Ranges
```python
# Check if data is in expected ranges
print("\nValue ranges:")
for col in ['DeNovoCNN_prob', 'DNM', 'MQ', 'VAF', 'child_coverage']:
    if col in ukbb_df.columns:
        print(f"  {col}: {ukbb_df[col].min():.2f} to {ukbb_df[col].max():.2f}")
```

### 3. Variant Type Column
```python
# Ensure var_type column exists
if 'var_type' not in ukbb_df.columns:
    print("⚠️  Missing 'var_type' column - needed for variant type filtering")
    # Should be created during data loading with values: 'SNV', 'Insertion', 'Deletion'
```

### 4. Run Diagnostic
```bash
python debug_optimization.py
```

This shows:
- Configured columns and their optimization types
- Expected parameter names
- Parameter name mappings

## Example: Complete Workflow

```python
from dnm_harmoniser import PipelineConfig, run_optimisation

# 1. Load configuration
config = PipelineConfig.from_yaml('filter.yaml')

# 2. Run optimization
result = run_optimisation(
    data='ukbb_data.tsv',
    reference='decode_reference.tsv',
    config=config
)

# 3. Get best parameters (NEW format)
all_best_params = result.best_params
# {'SNV': {'max_DeNovoCNN_prob': 0.5, ...}, 'Insertion': {...}, ...}

# 4. Use updated plotting function
from updated_plot_function import plot_results
plot_results(ukbb_df, decode_df, all_best_params)
```

## Key Takeaways

1. **Parameter naming**: Uses exact column names from config, not shortened versions
2. **Minimum optimization**: Counter-intuitively uses `max_*` thresholds (filter out high values)
3. **Maximum optimization**: Uses `min_*` thresholds (filter out low values)
4. **Range optimization**: Provides both `min_*` and `max_*` thresholds
5. **Symmetric ranges**: No manual calculation needed (e.g., `100 - min_vaf`)
6. **Linked columns**: Get separate parameters but same values

## Files

- [debug_optimization.py](debug_optimization.py) - Diagnostic tool
- [updated_plot_function.py](updated_plot_function.py) - Corrected plotting function
- [filter.yaml](/Users/kartikchundru/dnms/ukb/filter.yaml) - Your configuration
