# Fix: 'OptimisationConfig' object has no attribute 'column_names'

## Problem

The error occurred because the pipeline code was trying to access `config.optimisation.column_names`, but the `OptimisationConfig` class doesn't have this attribute. Instead, it has a `columns` attribute which is a list of `ColumnConfig` objects.

## Root Cause

The `_suggest_params()` method in `pipeline.py` (line 360) was written for an older version of the configuration system that had a simple `column_names` list. The new system has:
- `columns`: List of `ColumnConfig` objects with detailed settings
- Each column has: `name`, `dtype`, `optimisation`, `range_constraint`, `linked_to`, `variant_types`

## Solution

### 1. Updated Parameter Suggestion Logic

**File**: `src/dnm_harmoniser/pipeline.py` (lines 355-433)

**Before**: Used heuristics based on column names
```python
# Old approach - guess based on name
if 'DP' in col or 'depth' in col.lower():
    # Assume it's depth
elif 'VAF' in col:
    # Assume it's variant allele frequency
```

**After**: Uses explicit column configuration
```python
# New approach - use configuration
opt_columns = self.config.optimisation.get_optimisation_columns()

for col_config in opt_columns:
    if col_config.optimisation == 'minimum':
        # Suggest maximum threshold (filter values above)
    elif col_config.optimisation == 'maximum':
        # Suggest minimum threshold (filter values below)
    elif col_config.optimisation == 'range':
        # Suggest both min and max within constraints
```

### 2. Optimization Type Handling

The method now properly handles the three optimization types:

#### Minimum Optimization
For columns like `DeNovoCNN_prob`, `DNM`, `MQ` (where lower is better):
- Suggests a **maximum threshold**
- Values **above** this threshold are filtered out
- Parameter name: `max_{column_name}`

#### Maximum Optimization
For columns like `child_coverage`, `father_coverage` (where higher is better):
- Suggests a **minimum threshold**
- Values **below** this threshold are filtered out
- Parameter name: `min_{column_name}`

#### Range Optimization
For columns like `VAF`, `IMF` (where values should be within a range):
- Suggests both **min and max thresholds**
- Uses `range_constraint` bounds if specified
- Parameter names: `min_{column_name}` and `max_{column_name}`

### 3. Linked Columns Support

Updated linked columns handling to use configuration:
```python
linked_groups = self.config.optimisation.get_linked_column_groups()
for group in linked_groups:
    # Apply the first column's threshold to all linked columns
    first_col = group[0].name
    for linked_col in group[1:]:
        params[f'min_{linked_col.name}'] = params[f'min_{first_col.name}']
```

## Example from Configuration

From `/Users/kartikchundru/dnms/ukb/filter.yaml`:

```yaml
columns:
  # Minimum optimization - lower is better
  - name: DeNovoCNN_prob
    dtype: float
    optimisation: minimum  # → suggests max_DeNovoCNN_prob

  # Maximum optimization - higher is better
  - name: child_coverage
    dtype: int
    optimisation: maximum  # → suggests min_child_coverage

  # Linked columns - share same threshold
  - name: father_coverage
    dtype: int
    optimisation: maximum
    linked_to: mother_coverage  # → min_father_coverage = min_mother_coverage

  - name: mother_coverage
    dtype: int
    optimisation: maximum
    linked_to: father_coverage

  # Range optimization with symmetric constraint
  - name: VAF
    dtype: float
    optimisation: range
    range_constraint:
      lower: 25
      scale: 100  # → suggests min_VAF and max_VAF within [25, 75]
```

## How It Works Now

### 1. Configuration Defines Columns
```yaml
optimisation:
  columns:
    - name: MQ
      dtype: int
      optimisation: minimum  # Lower mapping quality is better
```

### 2. Pipeline Reads Configuration
```python
opt_columns = config.optimisation.get_optimisation_columns()
# Returns: [ColumnConfig(name='MQ', optimisation='minimum', dtype='int'), ...]
```

### 3. Optuna Suggests Thresholds
```python
for col_config in opt_columns:
    if col_config.optimisation == 'minimum':
        # MQ minimum = lower quality is better
        # So we suggest a MAXIMUM threshold
        params['max_MQ'] = trial.suggest_int('max_MQ', min_value, max_value)
        # Variants with MQ > max_MQ are filtered out
```

### 4. Filters Are Applied
```python
# In apply_filters()
if 'max_MQ' in params:
    data = data[data['MQ'] <= params['max_MQ']]  # Keep low MQ values
```

## Benefits

1. **Type-Safe**: Uses configuration objects instead of guessing from names
2. **Explicit**: Clear what each column's optimization goal is
3. **Flexible**: Supports all optimization types (minimum, maximum, range)
4. **Linked Columns**: Properly handles shared thresholds
5. **Range Constraints**: Respects configured bounds for range optimization
6. **Variant-Specific**: Ready to support variant-type-specific columns

## Testing

```python
from dnm_harmoniser import PipelineConfig

config = PipelineConfig.from_yaml('/Users/kartikchundru/dnms/ukb/filter.yaml')

# Get optimization columns
opt_cols = config.optimisation.get_optimisation_columns()
# Returns 9 columns configured for optimization

# Get linked groups
linked = config.optimisation.get_linked_column_groups()
# Returns: [['father_coverage', 'mother_coverage']]
```

✅ Configuration loads successfully
✅ Optimization columns identified correctly
✅ Linked columns grouped properly
✅ Parameter suggestion logic updated
✅ All optimization types supported

## Files Modified

1. **src/dnm_harmoniser/pipeline.py** (lines 355-433)
   - Removed `column_names` reference
   - Updated to use `get_optimisation_columns()`
   - Rewrote parameter suggestion logic
   - Updated linked columns handling
