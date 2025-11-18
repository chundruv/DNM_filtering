# Configuration Changelog

## Latest Changes

### Required Metadata Column Names (BREAKING CHANGE)

**What changed:**
- Sample ID, paternal age, and maternal age are no longer configured as regular columns
- Instead, you now specify the column names in your data using dedicated fields

**Why this change:**
- Clearer separation between required metadata and optimization columns
- Makes it explicit which columns are used for regression analysis
- Easier to configure when your data uses different column names

**Migration Guide:**

**Before (Old Format):**
```yaml
optimization:
  columns:
    - name: sample_id
      dtype: str
      optimization: none

    - name: paternal_age
      dtype: int
      optimization: none

    - name: maternal_age
      dtype: int
      optimization: none

    - name: quality_score
      dtype: float
      optimization: minimize
    # ... other columns
```

**After (New Format):**
```yaml
optimization:
  # Specify your column names
  sample_id_column: sample_id        # or "SampleID", "ID", etc.
  paternal_age_column: paternal_age  # or "Father_Age", etc.
  maternal_age_column: maternal_age  # or "Mother_Age", etc.
  reference_column: ref              # or "REF", etc.

  columns:
    - name: quality_score
      dtype: float
      optimization: minimize
    # ... other columns
```

### New Fields in OptimizationConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `sample_id_column` | str | `"sample_id"` | Column name for sample IDs in your data |
| `paternal_age_column` | str | `"paternal_age"` | Column name for father's age in your data |
| `maternal_age_column` | str | `"maternal_age"` | Column name for mother's age in your data |
| `reference_column` | str | `"ref"` | Column name for reference allele in your data |

### New Helper Method

```python
# Get the required metadata column names
required_meta = opt_config.get_required_metadata_columns()
# Returns: {
#   'sample_id': 'sample_id',
#   'paternal_age': 'paternal_age',
#   'maternal_age': 'maternal_age',
#   'reference': 'ref'
# }
```

### Updated Behavior

- `get_metadata_columns()` now returns only **additional** metadata columns (those with `optimization="none"`)
- Required metadata (sample_id, ages, reference) are accessed via `get_required_metadata_columns()`

## Previous Features

All previous features remain unchanged:

✓ Column types: `minimize`, `maximize`, `range`, `none`
✓ Linked columns (shared thresholds)
✓ Variant-specific columns (SNV, Insertion, Deletion)
✓ Range constraints with validation
✓ Full Pydantic validation

## Examples

### Example 1: Standard Column Names

If your data uses standard names, just use the defaults:

```yaml
optimization:
  sample_id_column: sample_id
  paternal_age_column: paternal_age
  maternal_age_column: maternal_age
  reference_column: ref
```

### Example 2: Custom Column Names

If your data uses different names:

```yaml
optimization:
  sample_id_column: SampleID
  paternal_age_column: Father_Age
  maternal_age_column: Mother_Age
  reference_column: REF
```

### Example 3: Python API

```python
from dnm_harmoniser.config import OptimizationConfig, ColumnConfig

config = OptimizationConfig(
    # Your column names
    sample_id_column="SampleID",
    paternal_age_column="Father_Age",
    maternal_age_column="Mother_Age",
    reference_column="REF",

    # Optimization columns
    columns=[
        ColumnConfig(name="quality", dtype="float", optimization="minimize"),
        ColumnConfig(name="depth", dtype="int", optimization="maximize"),
    ]
)

# Access the column names
meta = config.get_required_metadata_columns()
print(meta['sample_id'])  # Output: "SampleID"
```

## Compatibility

- ✅ YAML configurations must be updated to use the new format
- ✅ Python API users must specify metadata column names when creating OptimizationConfig
- ✅ All other features remain backward compatible
