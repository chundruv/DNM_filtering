# Column Configuration Guide

This guide explains how to configure columns in the DNM filtering pipeline.

## Required Metadata Column Names

First, specify the column names in your data for required metadata fields:

```yaml
optimization:
  # Column names in your data
  sample_id_column: sample_id        # or "SampleID", "ID", etc.
  paternal_age_column: paternal_age  # or "Father_Age", "dad_age", etc.
  maternal_age_column: maternal_age  # or "Mother_Age", "mom_age", etc.
  reference_column: ref              # or "REF", "reference", etc.
```

These fields tell the pipeline which columns in your data contain:
- **sample_id_column**: Sample identifiers
- **paternal_age_column**: Father's age (used in regression)
- **maternal_age_column**: Mother's age (used in regression)
- **reference_column**: Reference allele

## Column Types

### 1. Additional Metadata Columns (`optimization="none"`)
If you have other metadata columns to preserve (beyond the required ones above), add them with `optimization="none"`.

**Examples:**
- Other identifiers or annotations
- Additional demographic information

```yaml
columns:
  - name: population
    dtype: str
    optimization: none

  - name: sequencing_batch
    dtype: str
    optimization: none
```

### 2. Optimization Columns

#### Minimize (`optimization="minimize"`)
Lower values are better. Used for quality metrics where lower is preferred.

```yaml
- name: quality_score
  dtype: float
  optimization: minimize
```

#### Maximize (`optimization="maximize"`)
Higher values are better. Used for metrics like coverage or depth.

```yaml
- name: coverage
  dtype: int
  optimization: maximize
```

#### Range (`optimization="range"`)
Keep values within a specified range. Requires `range_constraint` with `min` and `max`.

```yaml
- name: allele_balance
  dtype: float
  optimization: range
  range_constraint:
    min: 0.25
    max: 0.75
```

### 3. Linked Columns

Linked columns share the same threshold value. This is useful when you want to apply the same cutoff to related columns (e.g., father and mother depth).

**Requirements:**
- Both columns must have the same `optimization` type
- Links must be bidirectional (both columns must reference each other)
- Cannot be used with metadata columns (`optimization="none"`)

```yaml
columns:
  - name: depth_father
    dtype: int
    optimization: maximize
    linked_to: depth_mother

  - name: depth_mother
    dtype: int
    optimization: maximize
    linked_to: depth_father
```

### 4. Variant-Specific Columns

Columns can be restricted to specific variant types using the `variant_types` field. If omitted, the column applies to all variant types.

**Use cases:**
- `base_quality`: Only relevant for SNVs
- `indel_length`: Only relevant for insertions and deletions
- `homopolymer_length`: Only relevant for deletions

```yaml
columns:
  # SNV-only column
  - name: base_quality
    dtype: float
    optimization: minimize
    variant_types:
      - SNV

  # Indel-only column (insertions and deletions)
  - name: indel_length
    dtype: int
    optimization: range
    range_constraint:
      min: 1
      max: 50
    variant_types:
      - Insertion
      - Deletion

  # Deletion-only column
  - name: homopolymer_length
    dtype: int
    optimization: maximize
    variant_types:
      - Deletion
```

## Complete Example

```yaml
optimization:
  # Required metadata column names
  sample_id_column: sample_id
  paternal_age_column: paternal_age
  maternal_age_column: maternal_age
  reference_column: ref

  columns:
    # Regular optimization columns (apply to all variant types)
    - name: quality_score
      dtype: float
      optimization: minimize

    - name: coverage
      dtype: int
      optimization: maximize

    - name: allele_balance
      dtype: float
      optimization: range
      range_constraint:
        min: 0.25
        max: 0.75

    # Variant-specific columns
    - name: base_quality
      dtype: float
      optimization: minimize
      variant_types:
        - SNV

    - name: indel_length
      dtype: int
      optimization: range
      range_constraint:
        min: 1
        max: 50
      variant_types:
        - Insertion
        - Deletion

    # Linked columns (share the same threshold)
    - name: depth_father
      dtype: int
      optimization: maximize
      linked_to: depth_mother

    - name: depth_mother
      dtype: int
      optimization: maximize
      linked_to: depth_father
```

## Python API

### Creating Configurations

```python
from dnm_harmoniser.config import ColumnConfig, RangeConstraint, OptimizationConfig

# Regular optimization column
opt_col = ColumnConfig(
    name="quality_score",
    dtype="float",
    optimization="minimize"
)

# Range-constrained column
range_col = ColumnConfig(
    name="allele_balance",
    dtype="float",
    optimization="range",
    range_constraint=RangeConstraint(min=0.25, max=0.75)
)

# Linked columns
father_depth = ColumnConfig(
    name="depth_father",
    dtype="int",
    optimization="maximize",
    linked_to="depth_mother"
)

mother_depth = ColumnConfig(
    name="depth_mother",
    dtype="int",
    optimization="maximize",
    linked_to="depth_father"
)

# Variant-specific columns
snv_only = ColumnConfig(
    name="base_quality",
    dtype="float",
    optimization="minimize",
    variant_types=["SNV"]
)

indel_only = ColumnConfig(
    name="indel_length",
    dtype="int",
    optimization="range",
    range_constraint=RangeConstraint(min=1, max=50),
    variant_types=["Insertion", "Deletion"]
)

# Create optimization config
opt_config = OptimizationConfig(
    # Specify your column names
    sample_id_column="SampleID",
    paternal_age_column="Father_Age",
    maternal_age_column="Mother_Age",
    reference_column="REF",
    # Add optimization columns
    columns=[opt_col, range_col, father_depth, mother_depth, snv_only, indel_only]
)
```

### Helper Methods

```python
# Get required metadata column names
required_meta = opt_config.get_required_metadata_columns()
print(required_meta)
# Output: {'sample_id': 'SampleID', 'paternal_age': 'Father_Age',
#          'maternal_age': 'Mother_Age', 'reference': 'REF'}

# Get all additional metadata columns (optimization='none')
additional_metadata = opt_config.get_metadata_columns()

# Get all columns that will be optimized
to_optimize = opt_config.get_optimization_columns()

# Get all columns for a specific variant type
snv_columns = opt_config.get_columns_for_variant_type("SNV")
insertion_columns = opt_config.get_columns_for_variant_type("Insertion")
deletion_columns = opt_config.get_columns_for_variant_type("Deletion")

# Get optimization columns (excluding metadata) for a specific variant type
snv_opt_columns = opt_config.get_optimization_columns_for_variant_type("SNV")

# Get groups of linked columns
linked_groups = opt_config.get_linked_column_groups()
for group in linked_groups:
    print(f"Linked: {[col.name for col in group]}")
```

## Validation

The configuration system automatically validates:

1. **Range constraints**: `min` must be less than `max`
2. **Range optimization**: `range_constraint` is required when `optimization="range"`
3. **Linked columns**:
   - Both columns must exist in the configuration
   - Both must have the same optimization type
   - Links must be bidirectional
4. **Metadata columns**: Cannot have `linked_to` set
5. **Variant types**: If specified, must contain at least one variant type

Any validation errors will be raised when creating or loading the configuration.
