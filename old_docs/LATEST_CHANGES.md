# Latest Configuration Changes

## Summary of Changes

### 1. Separate Data and Reference Column Names

**What changed:** All metadata columns now require two values - one for input data and one for reference data.

**Fields Updated:**
- `sample_id_column` → `sample_id_column_data` + `sample_id_column_reference`
- `paternal_age_column` → `paternal_age_column_data` + `paternal_age_column_reference`
- `maternal_age_column` → `maternal_age_column_data` + `maternal_age_column_reference`
- `reference_column` → `reference_column_data` + `reference_column_reference`

**New Field:**
- `alternate_column_data` + `alternate_column_reference` (for alternate alleles)

### 2. UK English Spelling

**What changed:** All American English spelling changed to British English.

**Spelling Changes:**
- `optimization` → `optimisation`
- `optimize` → `optimise`
- `minimize` → `minimum`
- `maximize` → `maximum`
- `OptimizationConfig` → `OptimisationConfig`

## New Configuration Format

### YAML Example

```yaml
optimisation:  # Changed from "optimization"
  variant_types:
    - SNV
    - Insertion
    - Deletion

  # Input data column names
  sample_id_column_data: sample_id
  paternal_age_column_data: paternal_age
  maternal_age_column_data: maternal_age
  reference_column_data: ref
  alternate_column_data: alt

  # Reference data column names
  sample_id_column_reference: sample_id
  paternal_age_column_reference: paternal_age
  maternal_age_column_reference: maternal_age
  reference_column_reference: ref
  alternate_column_reference: alt

  columns:
    - name: quality_score
      dtype: float
      optimisation: minimum  # Changed from "optimization: minimize"

    - name: coverage
      dtype: int
      optimisation: maximum  # Changed from "optimization: maximize"

    - name: allele_balance
      dtype: float
      optimisation: range
      range_constraint:
        min: 0.25
        max: 0.75
```

### Python API Example

```python
from dnm_harmoniser.config import OptimisationConfig, ColumnConfig  # Changed class name

config = OptimisationConfig(
    # Input data columns
    sample_id_column_data="SampleID",
    paternal_age_column_data="Father_Age",
    maternal_age_column_data="Mother_Age",
    reference_column_data="REF",
    alternate_column_data="ALT",

    # Reference data columns
    sample_id_column_reference="sample_id",
    paternal_age_column_reference="paternal_age",
    maternal_age_column_reference="maternal_age",
    reference_column_reference="ref",
    alternate_column_reference="alt",

    columns=[
        ColumnConfig(name="quality", dtype="float", optimisation="minimum"),  # UK spelling
        ColumnConfig(name="depth", dtype="int", optimisation="maximum"),      # UK spelling
    ]
)

# New helper methods
data_cols = config.get_required_metadata_columns_data()
ref_cols = config.get_required_metadata_columns_reference()
```

## New Helper Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `get_required_metadata_columns_data()` | `Dict[str, str]` | Column names for input data |
| `get_required_metadata_columns_reference()` | `Dict[str, str]` | Column names for reference data |
| `get_metadata_columns()` | `List[ColumnConfig]` | Additional metadata columns |
| `get_optimisation_columns()` | `List[ColumnConfig]` | Columns to be optimised |

## Migration Guide

### Step 1: Update Field Names

**Old:**
```yaml
optimisation:
  sample_id_column: sample_id
  paternal_age_column: paternal_age
  maternal_age_column: maternal_age
  reference_column: ref
```

**New:**
```yaml
optimisation:
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
```

### Step 2: Update Spelling

**Old:**
```yaml
columns:
  - name: quality_score
    dtype: float
    optimization: minimize
```

**New:**
```yaml
columns:
  - name: quality_score
    dtype: float
    optimisation: minimum
```

### Step 3: Update Python Imports

**Old:**
```python
from dnm_harmoniser.config import OptimizationConfig

config.get_required_metadata_columns()
```

**New:**
```python
from dnm_harmoniser.config import OptimisationConfig

data_cols = config.get_required_metadata_columns_data()
ref_cols = config.get_required_metadata_columns_reference()
```

## Key Benefits

1. **Separate Data Sources**: Clear distinction between input data and reference data column names
2. **Alternate Allele Support**: Explicit field for alternate allele columns
3. **UK English**: Consistent British English spelling throughout
4. **Better Flexibility**: Can use different column names in data vs reference

## Compatibility

⚠️ **BREAKING CHANGES** - Old configurations will need to be updated:

1. All `optimization` → `optimisation`
2. All `minimize`/`maximize` → `minimum`/`maximum`
3. Single column fields split into `_data` and `_reference` versions
4. New `alternate_column_data` and `alternate_column_reference` fields required

## Example: Different Column Names

If your input data uses different names than your reference:

```yaml
optimisation:
  # Input data uses uppercase
  sample_id_column_data: SampleID
  paternal_age_column_data: Father_Age
  maternal_age_column_data: Mother_Age
  reference_column_data: REF
  alternate_column_data: ALT

  # Reference uses lowercase
  sample_id_column_reference: sample_id
  paternal_age_column_reference: paternal_age
  maternal_age_column_reference: maternal_age
  reference_column_reference: ref
  alternate_column_reference: alt
```

## Testing

All changes have been tested and validated:
- ✅ YAML configuration loading
- ✅ Separate data/reference column names
- ✅ Alternate allele column support
- ✅ UK English spelling (optimisation, minimum, maximum)
- ✅ Helper methods working correctly
- ✅ All existing features preserved
