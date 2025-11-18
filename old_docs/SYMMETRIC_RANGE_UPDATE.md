# Symmetric Range Constraint - Implementation Summary

## What Was Implemented

A new `SymmetricRangeConstraint` class that allows you to specify symmetric range constraints where the upper bound is automatically calculated as `upper = scale - lower`.

## Key Features

1. **Automatic Upper Bound Calculation**: Specify only the lower bound, and the upper bound is calculated automatically
2. **Flexible Scaling**: Supports different scales (1.0 for 0-1 range, 100 for percentage)
3. **Full Validation**: Ensures the constraint makes mathematical sense
4. **Backward Compatible**: Regular `RangeConstraint` continues to work for non-symmetric ranges
5. **Easy to Use**: Works seamlessly in both YAML and Python API

## Use Cases

Perfect for metrics that are symmetric around a midpoint:
- **Allele Balance**: Keep between 0.25-0.75 (symmetric around 0.5)
- **VAF (Variant Allele Frequency)**: Keep between 25%-75% (symmetric around 50%)
- **Heterozygosity**: Any metric where you want to exclude extreme values on both sides

## Quick Examples

### YAML Configuration

**Allele Balance (0-1 scale):**
```yaml
- name: allele_balance
  dtype: float
  optimisation: range
  range_constraint:
    lower: 0.25
    scale: 1.0  # Results in upper = 0.75
```

**VAF Percentage (0-100 scale):**
```yaml
- name: vaf_percent
  dtype: float
  optimisation: range
  range_constraint:
    lower: 25
    scale: 100  # Results in upper = 75
```

### Python API

```python
from dnm_harmoniser.config import SymmetricRangeConstraint, ColumnConfig

# Create symmetric constraint
constraint = SymmetricRangeConstraint(lower=0.25, scale=1.0)
print(f"Range: [{constraint.min}, {constraint.max}]")  # [0.25, 0.75]

# Use in column configuration
col = ColumnConfig(
    name="allele_balance",
    dtype="float",
    optimisation="range",
    range_constraint=constraint
)
```

## Mathematical Relationship

```
upper = scale - lower

Examples:
- lower=0.25, scale=1.0  → upper=0.75
- lower=25,   scale=100  → upper=75
- lower=0.20, scale=1.0  → upper=0.80
- lower=30,   scale=100  → upper=70
```

## Validation

The constraint automatically validates:
1. `lower < scale/2` (lower must be less than midpoint)
2. `lower >= 0` (lower must be non-negative)
3. `upper <= scale` (upper doesn't exceed scale - automatically true if rule 1 passes)

## Files Changed

1. **[src/dnm_harmoniser/config.py](src/dnm_harmoniser/config.py)**
   - Added `SymmetricRangeConstraint` class (lines 130-175)
   - Updated `ColumnConfig.range_constraint` to accept `Union[RangeConstraint, SymmetricRangeConstraint]`

2. **[src/dnm_harmoniser/__init__.py](src/dnm_harmoniser/__init__.py)**
   - Added `SymmetricRangeConstraint` to exports

3. **[example_config.yaml](example_config.yaml)**
   - Updated `allele_balance` to use symmetric constraint
   - Added `vaf_percent` example with percentage scale

4. **Documentation**
   - Created [SYMMETRIC_RANGE_GUIDE.md](SYMMETRIC_RANGE_GUIDE.md) - Comprehensive guide
   - Updated [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - Added symmetric range examples

## Comparison: Before and After

### Before (Manual Upper Bound)

```yaml
- name: allele_balance
  dtype: float
  optimisation: range
  range_constraint:
    min: 0.25
    max: 0.75  # Must manually ensure 0.75 = 1.0 - 0.25
```

**Problems:**
- Easy to make mistakes (e.g., setting max: 0.70 instead of 0.75)
- Not clear that the range is intended to be symmetric
- Must manually calculate upper bound

### After (Automatic Upper Bound)

```yaml
- name: allele_balance
  dtype: float
  optimisation: range
  range_constraint:
    lower: 0.25
    scale: 1.0  # upper automatically calculated as 0.75
```

**Benefits:**
- Impossible to set non-symmetric bounds by mistake
- Clear semantic meaning (symmetric around midpoint)
- Less configuration to write
- Validated at load time

## Properties

The `SymmetricRangeConstraint` provides both its own properties and compatibility aliases:

| Property | Description |
|----------|-------------|
| `lower` | The lower bound (user-specified) |
| `upper` | The upper bound (calculated as `scale - lower`) |
| `scale` | The scale factor (default: 1.0) |
| `min` | Alias for `lower` (for compatibility with `RangeConstraint`) |
| `max` | Alias for `upper` (for compatibility with `RangeConstraint`) |

This means you can use `constraint.min` and `constraint.max` regardless of whether it's a `RangeConstraint` or `SymmetricRangeConstraint`.

## Testing

All functionality has been tested:
- ✅ Creating symmetric constraints (0-1 and 0-100 scale)
- ✅ Properties and aliases work correctly
- ✅ Using in ColumnConfig
- ✅ Loading from YAML
- ✅ Mixing symmetric and regular constraints
- ✅ Validation catches invalid constraints
- ✅ Backward compatibility with regular RangeConstraint

See [test_symmetric_constraint.py](test_symmetric_constraint.py) for full test suite.

## When to Use Each Type

### Use SymmetricRangeConstraint When:
- Range is symmetric around midpoint (e.g., 0.25-0.75 around 0.5)
- You want to exclude extreme values equally on both sides
- Examples: allele balance, VAF, heterozygosity

### Use Regular RangeConstraint When:
- Range is NOT symmetric (e.g., 40-60 is symmetric, but 35-65 around 50 is also valid)
- You need different distances from midpoint
- Examples: mapping quality, depth ranges, GC content

## Further Reading

- **Comprehensive Guide**: [SYMMETRIC_RANGE_GUIDE.md](SYMMETRIC_RANGE_GUIDE.md)
- **Quick Reference**: [QUICK_REFERENCE.md](QUICK_REFERENCE.md)
- **Column Configuration**: [COLUMN_CONFIG_GUIDE.md](COLUMN_CONFIG_GUIDE.md)
- **Example Config**: [example_config.yaml](example_config.yaml)
