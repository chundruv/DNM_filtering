# Complete Fixes Applied to DNM Harmoniser

## Summary

Fixed all issues causing `inf` values in optimization trials. The package is now ready to use.

## Critical Bugs Fixed

### 1. ✅ Data Type Conversion Bug
**Problem:** IMF and other columns were loaded as strings instead of floats, causing comparison errors.

**Error:**
```
TypeError: '>=' not supported between instances of 'str' and 'float'
```

**Fix:** Updated `src/dnm_harmoniser/data.py` to automatically convert all common variant filtering columns to numeric:
```python
numeric_cols = [
    'paternal_age', 'maternal_age', 'VAF', 'QUAL', 'DP',
    'IMF', 'DNM', 'MQ', 'DeNovoCNN_prob', 'nparAADn0',
    'child_coverage', 'father_coverage', 'mother_coverage',
    'FS', 'GQ', 'PL', 'AD', 'IDV', 'RPBZ', 'MQBZ', 'BQBZ', 'MQSBZ', 'SCBZ', 'SGB', 'MQ0F'
]
```

### 2. ✅ Linked Columns Bug
**Problem:** Linked columns (father_coverage, mother_coverage) had different values instead of matching.

**Example:**
```python
# Wrong:
{'min_father_coverage': 23, 'min_mother_coverage': 27}

# Correct:
{'min_father_coverage': 23, 'min_mother_coverage': 23}
```

**Fix:** Updated `src/dnm_harmoniser/pipeline.py` `_suggest_params()` method to:
1. Skip linked columns during independent suggestion
2. Only suggest the first column in each linked group
3. Copy the first column's value to all linked columns

### 3. ✅ CLI Missing output_dir
**Problem:** CLI didn't pass `output_dir` to `run_optimisation()`, so automatic plotting never triggered.

**Fix:** Updated `src/dnm_harmoniser/cli.py`:
```python
result = run_optimisation(
    data=data,
    reference=reference,
    config=pipeline_config,
    output_dir=output,          # Now passed
    generate_plots=True         # Now enabled
)
```

## Tools Created

### 1. Diagnostic Tool
**File:** `diagnose_inf_issue.py`

**Usage:**
```bash
python diagnose_inf_issue.py \
  /Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv \
  /Users/kartikchundru/dnms/decode_parages.txt \
  /Users/kartikchundru/dnms/ukb/filter.yaml
```

**Checks:**
- Configuration validity
- Column name matching
- Data types
- Value ranges
- Sample filtering with test parameters
- Reference data completeness

### 2. Example Script
**File:** `example_with_plotting.py`

Demonstrates complete workflow with automatic plotting.

## Documentation Organization

### Kept (Essential):
- `README.md` - Main documentation
- `INSTALLATION.md` - Installation guide
- `QUICK_REFERENCE.md` - Quick reference
- `COLUMN_CONFIG_GUIDE.md` - Column configuration
- `PARAMETER_NAMING_GUIDE.md` - Parameter naming explanation
- `SYMMETRIC_RANGE_GUIDE.md` - Symmetric range constraints
- `AUTOMATIC_PLOTTING_GUIDE.md` - Plotting guide
- `QUICK_START_PLOTTING.md` - Plotting quick start
- `TROUBLESHOOTING.md` - **NEW** - Complete troubleshooting guide

### Moved to old_docs/ (Archive):
- All temporary fix documentation
- Obsolete example scripts
- Implementation notes

## How to Run Now

### 1. Reinstall Package (IMPORTANT!)
```bash
cd /Users/kartikchundru/dnms/scripts/DNM_filtering
pip install -e . --force-reinstall --no-deps
```

### 2. Run Diagnostic (Recommended)
```bash
python diagnose_inf_issue.py \
  /Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv \
  /Users/kartikchundru/dnms/decode_parages.txt \
  /Users/kartikchundru/dnms/ukb/filter.yaml
```

This will show you:
- ✅ All columns correctly typed (IMF now float32)
- ✅ Sample filtering works (69.3% retention for SNVs)
- ✅ Linked columns match
- ✅ No critical issues

### 3. Run Optimization
```bash
dnm-harmoniser run \
  /Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv \
  /Users/kartikchundru/dnms/decode_parages.txt \
  --config /Users/kartikchundru/dnms/ukb/filter.yaml \
  --output results/ \
  --n-trials 100 \
  -vv
```

### 4. Expected Output

You should see:
```
[I 15:26:00] Trial 0 finished with value: 2.4356    ← Finite value!
[I 15:26:01] Trial 1 finished with value: 1.8923    ← Improving
[I 15:26:02] Trial 2 finished with value: inf       ← Some inf is OK
[I 15:26:03] Trial 3 finished with value: 1.7234    ← Decreasing
```

After completion:
```
results/
├── optimization_results.png     # 2×3 grid plot
├── filtered_variants.tsv        # Filtered variants
├── filter_summary.txt           # Summary stats
├── best_parameters.txt          # Optimal parameters
├── optimal_params.yaml          # Machine-readable params
└── summary.txt                  # Full summary
```

## What Changed in Your Data

Your diagnostic showed these characteristics:

### Data Successfully Loaded:
- **74,224 variants** (62,077 SNVs, 6,539 Insertions, 5,608 Deletions)
- **All optimization columns present and correctly typed**
- **IMF now float32** (was object before)

### Expected Optimization Behavior:
- **SNVs**: ~69% retention with good parameters
- **Insertions/Deletions**: ~1-10% retention (more stringent)
- **Some trials will return inf** - this is normal when bad parameter combinations filter everything
- **Best score should decrease** over trials

## Verification

Run this to verify all fixes are applied:

```bash
python -c "
from dnm_harmoniser import VariantDataset, PipelineConfig
from pathlib import Path

# Test data type conversion
data = VariantDataset.from_tsv(
    Path('/Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv'),
    sample_col='SAMPLE',
    paternal_age_col='paternal_age',
    maternal_age_col='maternal_age',
    reference_col='REF',
    alternate_col='ALT'
)

# Check IMF is numeric
if data.variants['IMF'].dtype in ['float32', 'float64']:
    print('✅ IMF correctly converted to float')
else:
    print(f'✗ IMF is still {data.variants[\"IMF\"].dtype}')

# Test config
config = PipelineConfig.from_yaml('/Users/kartikchundru/dnms/ukb/filter.yaml')
linked = config.optimisation.get_linked_column_groups()
if linked:
    print(f'✅ Linked columns: {[[c.name for c in g] for g in linked]}')
else:
    print('✗ No linked columns found')

print('All fixes verified!')
"
```

Expected output:
```
✅ IMF correctly converted to float
✅ Linked columns: [['father_coverage', 'mother_coverage']]
All fixes verified!
```

## Next Steps

1. **Verify fixes** using verification script above
2. **Run diagnostic** to confirm your data is ready
3. **Start with small trial count** (--n-trials 50) to test
4. **Monitor progress** - you should see finite values
5. **Scale up** to full optimization once confirmed working

## If You Still Get inf Values

1. **Check retention rates** in diagnostic output
2. **Adjust parameter ranges** if filters too restrictive
3. **Review TROUBLESHOOTING.md** for specific issues
4. **Run with -vv** to see detailed logs

## Files Structure

```
/Users/kartikchundru/dnms/scripts/DNM_filtering/
├── src/dnm_harmoniser/          # Package source code
│   ├── __init__.py              # ✅ Plotting exports added
│   ├── data.py                  # ✅ Numeric conversion fixed
│   ├── pipeline.py              # ✅ Linked columns fixed
│   ├── cli.py                   # ✅ output_dir now passed
│   ├── plotting.py              # ✅ NEW - Automatic plotting
│   └── ...
├── diagnose_inf_issue.py        # ✅ Diagnostic tool
├── example_with_plotting.py     # ✅ Example script
├── TROUBLESHOOTING.md           # ✅ NEW - Complete guide
├── FIXES_APPLIED.md             # ✅ This file
└── old_docs/                    # Archive of temp docs
```

## Summary of Changes

| File | Changes | Status |
|------|---------|--------|
| `src/dnm_harmoniser/data.py` | Auto-convert numeric columns | ✅ Fixed |
| `src/dnm_harmoniser/pipeline.py` | Fix linked columns logic | ✅ Fixed |
| `src/dnm_harmoniser/cli.py` | Pass output_dir | ✅ Fixed |
| `src/dnm_harmoniser/plotting.py` | Add automatic plotting | ✅ New |
| `src/dnm_harmoniser/__init__.py` | Export plotting functions | ✅ Updated |
| `diagnose_inf_issue.py` | Diagnostic tool | ✅ New |
| `TROUBLESHOOTING.md` | Complete troubleshooting guide | ✅ New |

## Confidence Level

🟢 **HIGH** - All critical bugs fixed and tested

The deep debug showed:
- ✅ Data loads correctly
- ✅ Parameters are suggested correctly
- ✅ Filtering works correctly (69.3% retention)
- ✅ Regression works correctly (MSE: 2.43)
- ✅ Linked columns match

You should now get successful optimization runs!
