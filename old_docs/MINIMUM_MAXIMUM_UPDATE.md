# Update: Changed "minimise/maximise" to "minimum/maximum"

## Summary

Changed all instances of "minimise" and "maximise" to "minimum" and "maximum" to better reflect what the optimisation strategies are doing.

## Files Updated

### 1. ✅ src/dnm_harmoniser/config.py
- **Literal type definition** (line 134): `Literal["minimum", "maximum", "range", "none"]`
- **Docstring examples**: Updated all example code
- **PRESETS**: Updated all preset configurations

### 2. ✅ example_config.yaml
- All `optimisation: minimise` → `optimisation: minimum`
- All `optimisation: maximise` → `optimisation: maximum`
- Comments updated

### 3. ✅ Documentation Files
- **LATEST_CHANGES.md**: All instances updated
- **COLUMN_CONFIG_GUIDE.md**: All instances updated
- **QUICK_REFERENCE.md**: All instances updated
- **CHANGELOG.md**: All instances updated

## New Values

### Optimisation Types

| Value | Description | Use Case |
|-------|-------------|----------|
| `minimum` | Lower is better | Quality scores, error rates |
| `maximum` | Higher is better | Coverage, depth, mapping quality |
| `range` | Keep within bounds | Allele balance, GC content |
| `none` | No thresholding | Metadata columns |

## Example Usage

### YAML Configuration
```yaml
columns:
  - name: quality_score
    dtype: float
    optimisation: minimum  # Changed from "minimise"

  - name: coverage
    dtype: int
    optimisation: maximum  # Changed from "maximise"

  - name: allele_balance
    dtype: float
    optimisation: range
    range_constraint:
      min: 0.25
      max: 0.75
```

### Python API
```python
from dnm_harmoniser.config import ColumnConfig

# Minimum
min_col = ColumnConfig(
    name="quality_score",
    dtype="float",
    optimisation="minimum"  # Changed from "minimise"
)

# Maximum
max_col = ColumnConfig(
    name="coverage",
    dtype="int",
    optimisation="maximum"  # Changed from "maximise"
)
```

## Testing

All changes have been tested and validated:
- ✅ YAML configuration loads correctly with minimum/maximum
- ✅ Python API accepts minimum/maximum values
- ✅ Type checking passes (Literal type updated)
- ✅ All preset configurations work
- ✅ Documentation updated consistently

## Search Results

Confirmed no remaining instances of "minimise" or "maximise" in:
- src/dnm_harmoniser/config.py
- example_config.yaml
- Documentation files

All occurrences have been successfully updated to "minimum" and "maximum".
