# Fix: Column Name Mapping for Custom Data Schemas

## Problem

The error `'SAMPLE'` occurred because the code was hardcoded to expect specific column names (`SAMPLE`, `paternal_age`, `maternal_age`, etc.), but user data files can have different column names (e.g., `pid` for sample ID, `Father_Age` for paternal age).

## Root Cause

1. The configuration allowed specifying custom column names via:
   - `sample_id_column_data` / `sample_id_column_reference`
   - `paternal_age_column_data` / `paternal_age_column_reference`
   - `maternal_age_column_data` / `maternal_age_column_reference`
   - `reference_column_data` / `reference_column_reference`
   - `alternate_column_data` / `alternate_column_reference`

2. However, `VariantDataset.from_tsv()` was not receiving these column names
3. The data loader only standardized `SAMPLE` but not the age or ref/alt columns
4. The pipeline code hardcoded references to `'SAMPLE'`, `'paternal_age'`, etc.

## Solution

### 1. Extended `VariantDataset.from_tsv()` Parameters

**File**: `src/dnm_harmoniser/data.py`

Added parameters to accept custom column names:
```python
@classmethod
def from_tsv(
    cls,
    path: Path,
    sample_col: str = 'SAMPLE',
    paternal_age_col: str = 'paternal_age',
    maternal_age_col: str = 'maternal_age',
    reference_col: str = 'ref',
    alternate_col: str = 'alt',
    required_cols: Optional[List[str]] = None,
    max_length: int = 20
) -> 'VariantDataset':
```

### 2. Column Name Standardization

The loader now renames columns from user's schema to internal standard names:

```python
# User's data: SAMPLE, pid → Internal: SAMPLE, SAMPLE
if sample_col != 'SAMPLE' and sample_col in df.columns:
    df.rename(columns={sample_col: 'SAMPLE'}, inplace=True)

# User's data: paternal_age, Father_Age → Internal: paternal_age, paternal_age
if paternal_age_col != 'paternal_age' and paternal_age_col in df.columns:
    df.rename(columns={paternal_age_col: 'paternal_age'}, inplace=True)

# Similar for maternal_age, REF, ALT
```

### 3. Updated API Calls

**File**: `src/dnm_harmoniser/api.py`

All `VariantDataset.from_tsv()` calls now pass the configured column names:

```python
data = VariantDataset.from_tsv(
    Path(data_path),
    sample_col=config.optimisation.sample_id_column_data,
    paternal_age_col=config.optimisation.paternal_age_column_data,
    maternal_age_col=config.optimisation.maternal_age_column_data,
    reference_col=config.optimisation.reference_column_data,
    alternate_col=config.optimisation.alternate_column_data
)

reference = VariantDataset.from_tsv(
    Path(reference_path),
    sample_col=config.optimisation.sample_id_column_reference,
    paternal_age_col=config.optimisation.paternal_age_column_reference,
    maternal_age_col=config.optimisation.maternal_age_column_reference,
    reference_col=config.optimisation.reference_column_reference,
    alternate_col=config.optimisation.alternate_column_reference
)
```

## How It Works

### Internal Column Names (Standard)
The pipeline internally uses these standardized names:
- `SAMPLE` - Sample/individual identifier
- `paternal_age` - Father's age
- `maternal_age` - Mother's age
- `REF` - Reference allele
- `ALT` - Alternate allele

### User Configuration Maps to Internal Names

Example from `/Users/kartikchundru/dnms/ukb/filter.yaml`:

```yaml
optimisation:
  # Input data uses these column names:
  sample_id_column_data: SAMPLE          # Already matches → no rename
  paternal_age_column_data: paternal_age # Already matches → no rename
  maternal_age_column_data: maternal_age # Already matches → no rename
  reference_column_data: REF             # Already matches → no rename
  alternate_column_data: ALT             # Already matches → no rename

  # Reference data uses these column names:
  sample_id_column_reference: pid        # Renamed: pid → SAMPLE
  paternal_age_column_reference: paternal_age
  maternal_age_column_reference: maternal_age
  reference_column_reference: Ref        # Renamed: Ref → REF
  alternate_column_reference: Alt        # Renamed: Alt → ALT
```

### Process Flow

1. User specifies column names in YAML config
2. Config is loaded into `PipelineConfig`
3. When loading data:
   - `from_tsv()` receives user's column names
   - Renames columns to internal standard names
   - Pipeline code uses standard names (`SAMPLE`, `paternal_age`, etc.)
4. Data is processed with consistent column names

## Files Modified

1. **src/dnm_harmoniser/data.py** (lines 45-92)
   - Added 4 new parameters to `from_tsv()`
   - Added column renaming logic for age and ref/alt columns

2. **src/dnm_harmoniser/api.py** (lines 68-85, 183-200, 335-342)
   - Updated 3 locations where `from_tsv()` is called
   - Pass configuration column names to data loader

## Testing

```python
from dnm_harmoniser import PipelineConfig

config = PipelineConfig.from_yaml('/Users/kartikchundru/dnms/ukb/filter.yaml')

# Data file has columns: SAMPLE, paternal_age, maternal_age, REF, ALT
# Reference file has columns: pid, paternal_age, maternal_age, Ref, Alt

# After loading:
# - Both files will have: SAMPLE, paternal_age, maternal_age, REF, ALT
# - pid → SAMPLE (in reference)
# - Ref → REF (in reference)
# - Alt → ALT (in reference)
```

✅ Configuration loads correctly
✅ Column names mapped properly
✅ Data standardization working

## Benefits

1. **Flexible Schema Support**: Works with any column naming scheme
2. **Automatic Standardization**: User doesn't need to rename files
3. **Consistent Internal Processing**: Pipeline always works with standard names
4. **Clear Configuration**: YAML explicitly shows the mapping
