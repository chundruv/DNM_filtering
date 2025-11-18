# Quick Reference Card

## Installation

```bash
pip install -e .
```

## Import

```python
from dnm_harmoniser.config import (
    PipelineConfig,
    ColumnConfig,
    RangeConstraint,
    SymmetricRangeConstraint,
    OptimisationConfig
)
```

## Required Metadata Column Names

First, specify the column names in your data:

| Field | Default | Description |
|-------|---------|-------------|
| `sample_id_column` | `"sample_id"` | Column name for sample IDs |
| `paternal_age_column` | `"paternal_age"` | Column name for father's age |
| `maternal_age_column` | `"maternal_age"` | Column name for mother's age |
| `reference_column` | `"ref"` | Column name for reference allele |

## Column Configuration Options

### Basic Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `name` | str | ✓ | Column name in dataframe |
| `dtype` | `"int"`, `"float"`, `"str"` | ✓ | Data type |
| `optimization` | `"minimize"`, `"maximize"`, `"range"`, `"none"` | ✓ | Optimization strategy |
| `range_constraint` | `RangeConstraint` | * | Min/max bounds (required if `optimization="range"`) |
| `linked_to` | str | ✗ | Name of column to link with |
| `variant_types` | List | ✗ | `["SNV"]`, `["Insertion"]`, `["Deletion"]`, or combinations |

### Optimization Types

| Type | Description | Use Case |
|------|-------------|----------|
| `minimize` | Lower is better | Quality scores, error rates |
| `maximize` | Higher is better | Coverage, depth, mapping quality |
| `range` | Keep within bounds | Allele balance, GC content |
| `none` | No thresholding | Metadata (sample_id, ages) |

### Variant Types

| Value | Description |
|-------|-------------|
| `None` (default) | Applies to all variant types |
| `["SNV"]` | Single nucleotide variants only |
| `["Insertion"]` | Insertions only |
| `["Deletion"]` | Deletions only |
| `["Insertion", "Deletion"]` | Both indel types |

## Common Patterns

### 1. Specifying Required Metadata Columns
```yaml
optimization:
  sample_id_column: SampleID
  paternal_age_column: Father_Age
  maternal_age_column: Mother_Age
  reference_column: REF
```

### 2. Additional Metadata Column
```python
ColumnConfig(name="population", dtype="str", optimization="none")
```

### 3. Minimize
```python
ColumnConfig(name="quality_score", dtype="float", optimization="minimize")
```

### 3. Maximize
```python
ColumnConfig(name="coverage", dtype="int", optimization="maximize")
```

### 4. Range (Regular)
```python
ColumnConfig(
    name="mapping_quality",
    dtype="int",
    optimisation="range",
    range_constraint=RangeConstraint(min=40, max=60)
)
```

### 4b. Range (Symmetric)
For symmetric ranges where upper = scale - lower:
```python
ColumnConfig(
    name="allele_balance",
    dtype="float",
    optimisation="range",
    range_constraint=SymmetricRangeConstraint(lower=0.25, scale=1.0)  # upper=0.75
)
```

### 5. Linked Columns
```python
# Must be bidirectional!
ColumnConfig(name="depth_father", dtype="int", optimization="maximize", linked_to="depth_mother")
ColumnConfig(name="depth_mother", dtype="int", optimization="maximize", linked_to="depth_father")
```

### 6. Variant-Specific
```python
# SNV only
ColumnConfig(name="base_quality", dtype="float", optimization="minimize", variant_types=["SNV"])

# Indels only
ColumnConfig(
    name="indel_length",
    dtype="int",
    optimization="range",
    range_constraint=RangeConstraint(min=1, max=50),
    variant_types=["Insertion", "Deletion"]
)
```

## Helper Methods

```python
config = PipelineConfig.from_yaml("config.yaml")
opt = config.optimization

# Get required metadata column names
required_metadata = opt.get_required_metadata_columns()
# Returns: {'sample_id': 'sample_id', 'paternal_age': 'paternal_age', ...}

# Get different column types
additional_metadata = opt.get_metadata_columns()  # Additional metadata columns
to_optimize = opt.get_optimization_columns()
linked_groups = opt.get_linked_column_groups()

# Get columns for specific variant type
snv_cols = opt.get_columns_for_variant_type("SNV")
snv_opt_cols = opt.get_optimization_columns_for_variant_type("SNV")
```

## CLI Commands

```bash
# Show help
dnm-harmoniser --help

# Initialize configuration
dnm-harmoniser init --output config.yaml

# Show presets
dnm-harmoniser show-presets

# Validate data
dnm-harmoniser validate data.tsv --config config.yaml

# Run optimization
dnm-harmoniser run data.tsv --config config.yaml --output results/
```

## YAML Template

```yaml
optimization:
  # Required metadata column names (specify your column names)
  sample_id_column: sample_id
  paternal_age_column: paternal_age
  maternal_age_column: maternal_age
  reference_column: ref

  columns:
    # Minimize
    - name: quality_score
      dtype: float
      optimization: minimize

    # Maximize
    - name: coverage
      dtype: int
      optimization: maximize

    # Symmetric range (upper = 1.0 - lower)
    - name: allele_balance
      dtype: float
      optimisation: range
      range_constraint:
        lower: 0.25
        scale: 1.0  # upper = 0.75

    # Regular range (non-symmetric)
    - name: mapping_quality
      dtype: int
      optimisation: range
      range_constraint:
        min: 40
        max: 60

    # Variant-specific
    - name: base_quality
      dtype: float
      optimization: minimize
      variant_types:
        - SNV

    # Linked columns
    - name: depth_father
      dtype: int
      optimization: maximize
      linked_to: depth_mother

    - name: depth_mother
      dtype: int
      optimization: maximize
      linked_to: depth_father
```

## Validation Rules

✓ `range_constraint.min` < `range_constraint.max`
✓ `range_constraint` required when `optimization="range"`
✓ Linked columns must exist and have same optimization type
✓ Links must be bidirectional
✓ Metadata columns (`optimization="none"`) cannot be linked
✓ `variant_types` must not be empty if specified

## Files

- **Configuration:** [example_config.yaml](example_config.yaml)
- **Full Guide:** [COLUMN_CONFIG_GUIDE.md](COLUMN_CONFIG_GUIDE.md)
- **Installation:** [INSTALLATION.md](INSTALLATION.md)
- **Code:** [src/dnm_harmoniser/config.py](src/dnm_harmoniser/config.py)
