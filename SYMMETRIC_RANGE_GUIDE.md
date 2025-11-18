# Symmetric Range Constraint Guide

## Overview

The `SymmetricRangeConstraint` allows you to specify range constraints where the upper and lower bounds are symmetric around a midpoint. This is particularly useful for metrics like **allele balance** or **VAF (Variant Allele Frequency)** where you want to filter values that are too extreme on either side.

## Mathematical Relationship

In a symmetric range constraint:
```
upper = scale - lower
```

Where:
- `lower`: The lower bound you specify
- `upper`: Automatically calculated upper bound
- `scale`: The maximum value of your data range
  - Use `1.0` for 0-1 normalized data (e.g., allele balance)
  - Use `100` for percentage data (e.g., VAF%)

## Examples

### Example 1: Allele Balance (0-1 scale)

If you want to keep allele balance between 0.25 and 0.75:

**YAML:**
```yaml
- name: allele_balance
  dtype: float
  optimisation: range
  range_constraint:
    lower: 0.25
    scale: 1.0  # Results in upper = 1.0 - 0.25 = 0.75
```

**Python:**
```python
from dnm_harmoniser.config import SymmetricRangeConstraint, ColumnConfig

constraint = SymmetricRangeConstraint(lower=0.25, scale=1.0)
print(constraint.upper)  # Output: 0.75

col = ColumnConfig(
    name="allele_balance",
    dtype="float",
    optimisation="range",
    range_constraint=constraint
)
```

### Example 2: VAF Percentage (0-100 scale)

If you want to keep VAF between 25% and 75%:

**YAML:**
```yaml
- name: vaf_percent
  dtype: float
  optimisation: range
  range_constraint:
    lower: 25
    scale: 100  # Results in upper = 100 - 25 = 75
```

**Python:**
```python
constraint = SymmetricRangeConstraint(lower=25, scale=100)
print(constraint.upper)  # Output: 75
```

### Example 3: Different Symmetric Points

You can specify different lower bounds to create different symmetric ranges:

**20-80 range (0-1 scale):**
```yaml
range_constraint:
  lower: 0.20
  scale: 1.0  # upper = 0.80
```

**30-70 range (percentage):**
```yaml
range_constraint:
  lower: 30
  scale: 100  # upper = 70
```

## Comparison with Regular RangeConstraint

### Regular RangeConstraint
For non-symmetric ranges where you specify both bounds explicitly:

**YAML:**
```yaml
- name: mapping_quality
  dtype: int
  optimisation: range
  range_constraint:
    min: 40
    max: 60  # Not symmetric - you specify both values
```

### SymmetricRangeConstraint
For symmetric ranges where upper is calculated from lower:

**YAML:**
```yaml
- name: allele_balance
  dtype: float
  optimisation: range
  range_constraint:
    lower: 0.25
    scale: 1.0  # upper calculated automatically as 0.75
```

## Properties and Methods

The `SymmetricRangeConstraint` class provides:

| Property | Type | Description |
|----------|------|-------------|
| `lower` | float/int | The lower bound (specified by user) |
| `upper` | float/int | The upper bound (calculated as `scale - lower`) |
| `scale` | float/int | The scale factor (default: 1.0) |
| `min` | float/int | Alias for `lower` (for compatibility) |
| `max` | float/int | Alias for `upper` (for compatibility) |

## Validation Rules

The symmetric range constraint validates:

1. **Lower < Midpoint**: The lower bound must be less than the midpoint (`scale / 2`)
   ```python
   # ✗ INVALID - lower >= midpoint
   SymmetricRangeConstraint(lower=0.6, scale=1.0)  # midpoint is 0.5

   # ✓ VALID
   SymmetricRangeConstraint(lower=0.4, scale=1.0)
   ```

2. **Non-negative Lower**: The lower bound must be >= 0
   ```python
   # ✗ INVALID
   SymmetricRangeConstraint(lower=-0.1, scale=1.0)

   # ✓ VALID
   SymmetricRangeConstraint(lower=0.1, scale=1.0)
   ```

3. **Upper <= Scale**: The calculated upper bound must not exceed the scale (automatically ensured by rule 1)

## Use Cases

### When to Use SymmetricRangeConstraint

Use `SymmetricRangeConstraint` when:
- Values should be symmetric around a midpoint (0.5 for 0-1 scale, 50 for percentage)
- You want to filter extreme values on both sides
- Common examples: allele balance, VAF, heterozygosity

### When to Use Regular RangeConstraint

Use regular `RangeConstraint` when:
- Bounds are not symmetric
- You need different distances from the midpoint
- Common examples: mapping quality (40-60), depth ranges, GC content

## Python API Examples

### Creating Symmetric Constraints

```python
from dnm_harmoniser.config import (
    SymmetricRangeConstraint,
    ColumnConfig,
    OptimisationConfig
)

# Create symmetric constraint
sym_constraint = SymmetricRangeConstraint(lower=0.25, scale=1.0)

# Use in column config
col = ColumnConfig(
    name="allele_balance",
    dtype="float",
    optimisation="range",
    range_constraint=sym_constraint
)

# Access bounds
print(f"Range: [{sym_constraint.min}, {sym_constraint.max}]")
# Output: Range: [0.25, 0.75]

# Use in optimization config
opt_config = OptimisationConfig(
    sample_id_column_data="sample_id",
    paternal_age_column_data="paternal_age",
    maternal_age_column_data="maternal_age",
    reference_column_data="ref",
    alternate_column_data="alt",
    sample_id_column_reference="sample_id",
    paternal_age_column_reference="paternal_age",
    maternal_age_column_reference="maternal_age",
    reference_column_reference="ref",
    alternate_column_reference="alt",
    columns=[col]
)
```

### Mixing Symmetric and Regular Constraints

```python
columns = [
    # Symmetric constraint for allele balance
    ColumnConfig(
        name="allele_balance",
        dtype="float",
        optimisation="range",
        range_constraint=SymmetricRangeConstraint(lower=0.25, scale=1.0)
    ),
    # Regular constraint for mapping quality
    ColumnConfig(
        name="mapping_quality",
        dtype="int",
        optimisation="range",
        range_constraint=RangeConstraint(min=40, max=60)
    )
]
```

## Complete YAML Example

```yaml
optimisation:
  variant_types:
    - SNV
    - Insertion
    - Deletion

  # Required metadata
  sample_id_column_data: sample_id
  paternal_age_column_data: paternal_age
  maternal_age_column_data: maternal_age
  reference_column_data: ref
  alternate_column_data: alt

  sample_id_column_reference: sample_id
  paternal_age_column_reference: paternal_age
  maternal_age_column_reference: maternal_age
  reference_column_reference: ref
  alternate_column_reference: alt

  columns:
    # Symmetric range - allele balance (0-1 scale)
    - name: allele_balance
      dtype: float
      optimisation: range
      range_constraint:
        lower: 0.25
        scale: 1.0  # upper = 0.75

    # Symmetric range - VAF percentage (0-100 scale)
    - name: vaf_percent
      dtype: float
      optimisation: range
      range_constraint:
        lower: 25
        scale: 100  # upper = 75

    # Regular (non-symmetric) range
    - name: mapping_quality
      dtype: int
      optimisation: range
      range_constraint:
        min: 40
        max: 60

    # Other optimization types
    - name: quality_score
      dtype: float
      optimisation: minimum

    - name: coverage
      dtype: int
      optimisation: maximum
```

## Benefits

1. **Automatic Calculation**: Upper bound is automatically calculated, reducing configuration errors
2. **Semantic Clarity**: Makes it clear that the range is symmetric around the midpoint
3. **Validation**: Ensures the constraint makes mathematical sense
4. **Flexibility**: Can use different scales (0-1, 0-100, etc.)
5. **Backward Compatible**: Regular `RangeConstraint` still works for non-symmetric ranges
