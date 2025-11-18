# Fixing the "inf" Trial Values Issue

## Root Cause Identified

The primary issue was that the `IMF` column (and potentially other columns) were being loaded as **strings** (`object` dtype) instead of **floats**. This caused the filtering comparisons to fail with:

```
TypeError: '>=' not supported between instances of 'str' and 'float'
```

## Fix Applied

Updated `src/dnm_harmoniser/data.py` to automatically convert all common variant filtering columns to numeric types:

```python
numeric_cols = [
    'paternal_age', 'maternal_age', 'VAF', 'QUAL', 'DP',
    'IMF', 'DNM', 'MQ', 'DeNovoCNN_prob', 'nparAADn0',
    'child_coverage', 'father_coverage', 'mother_coverage',
    'FS', 'GQ', 'PL', 'AD', 'IDV', 'RPBZ', 'MQBZ', 'BQBZ', 'MQSBZ', 'SCBZ', 'SGB', 'MQ0F'
]
```

All these columns are now automatically converted to numeric during data loading.

## Diagnostic Results

Running the diagnostic script on your data shows:

### ✓ Data Successfully Loaded

- **Input data**: 74,224 variants (62,077 SNVs, 6,539 Insertions, 5,608 Deletions)
- **Reference data**: 381,578 variants (349,110 SNVs, 10,462 Insertions, 22,006 Deletions)
- **All optimization columns present** and correctly typed

### ✓ Columns Correctly Typed

| Column | Type | Non-null | Range |
|--------|------|----------|-------|
| DeNovoCNN_prob | float32 | 74224/74224 | 0.901 - 1.000 |
| DNM | float32 | 74224/74224 | -0.335 - 0.000 |
| MQ | int8 | 74224/74224 | 3 - 60 |
| nparAADn0 | int8 | 74224/74224 | 0 - 2 |
| child_coverage | int16 | 74224/74224 | 0 - 243 |
| father_coverage | int16 | 74224/74224 | 0 - 2339 |
| mother_coverage | int16 | 74224/74224 | 1 - 1052 |
| VAF | int8 | 74224/74224 | 14 - 93 |
| **IMF** | **float32** | 12147/74224 | 0.095 - 1.000 |

Note: IMF has NaN values for SNVs (as expected - it's only for indels).

### Sample Filtering Test Results

| Variant Type | Original | Filtered | Retention |
|--------------|----------|----------|-----------|
| SNV | 62,077 | 1,966 | 3.2% |
| Insertion | 6,539 | 77 | 1.2% |
| Deletion | 5,608 | 54 | 1.0% |

**Note**: The low retention rates with median parameter values suggest that optimal parameters will need careful tuning by the optimizer.

## Running the Optimization

Now that the data type issue is fixed, you can run the optimization:

### Method 1: Using CLI (Recommended)

```bash
dnm-harmoniser run \
  /Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv \
  /Users/kartikchundru/dnms/decode_parages.txt \
  --config /Users/kartikchundru/dnms/ukb/filter.yaml \
  --output results/ \
  -vv
```

The `-vv` flag provides detailed logging to help diagnose any remaining issues.

### Method 2: Using Python API

```python
from dnm_harmoniser import run_optimisation
from pathlib import Path

result = run_optimisation(
    data="/Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv",
    reference="/Users/kartikchundru/dnms/decode_parages.txt",
    config_file="/Users/kartikchundru/dnms/ukb/filter.yaml",
    output_dir="results/",
    generate_plots=True
)

print(result.summary)
```

## Expected Output Files

After successful optimization, you'll find in `results/`:

```
results/
├── optimization_results.png     # 2×3 grid plot (paternal/maternal × SNV/Ins/Del)
├── filtered_variants.tsv        # Filtered variants
├── filter_summary.txt           # Summary statistics
├── best_parameters.txt          # Optimal parameters
├── optimal_params.yaml          # Machine-readable parameters
├── summary.txt                  # Optimization summary
└── config.yaml                  # Copy of configuration used
```

## If You Still Get "inf" Values

If optimization still returns `inf` for all trials, it means the parameter combinations are filtering out ALL variants. To debug:

### 1. Run with verbose logging

```bash
dnm-harmoniser run ... -vv > optimization.log 2>&1
```

Check the log for:
- Which parameters are being suggested
- How many variants remain after filtering
- Any error messages

### 2. Adjust optimization ranges

Your config file might have parameter ranges that are too restrictive. Consider:

**For "minimum" optimization columns** (lower is better):
- The optimizer suggests `max_*` thresholds
- If the range is too low, all variants will be filtered out

**For "maximum" optimization columns** (higher is better):
- The optimizer suggests `min_*` thresholds
- If the range is too high, all variants will be filtered out

**For "range" optimization columns**:
- Both min and max are suggested
- If the range is too narrow, all variants will be filtered out

### 3. Start with fewer trials

Use a smaller number of trials to test:

```bash
dnm-harmoniser run ... --n-trials 10 -vv
```

### 4. Check your data distribution

Run the diagnostic script to understand your data:

```bash
python diagnose_inf_issue.py \
  /Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv \
  /Users/kartikchundru/dnms/decode_parages.txt \
  /Users/kartikchundru/dnms/ukb/filter.yaml
```

This shows:
- Data ranges for each column
- Sample filtering with median parameters
- Any missing or incorrectly typed columns

## Key Parameter Ranges from Your Data

Based on the diagnostic, here are the actual ranges in your data:

| Column | Min | Median | Max | Note |
|--------|-----|--------|-----|------|
| DeNovoCNN_prob | 0.901 | 0.990 | 1.000 | Very high values - most variants are confident |
| DNM | -0.335 | 0.000 | 0.000 | Most values are exactly 0 |
| MQ | 3 | 60 | 60 | Most values are maximum (60) |
| nparAADn0 | 0 | 0 | 2 | Most values are 0 |
| child_coverage | 0 | 31 | 243 | Good spread |
| VAF | 14 | 48 | 93 | Good spread |
| IMF | 0.095 | 0.404 | 1.000 | Only for indels |

### Implications for Optimization

1. **DeNovoCNN_prob**: Since median is 0.990, the optimizer will likely suggest max thresholds close to 1.0
2. **DNM**: Since most values are 0, filtering may be very stringent
3. **MQ**: Since most values are 60 (maximum), this may not be a good discriminator
4. **nparAADn0**: Since most values are 0, filtering on this will be binary

## Monitoring Optimization Progress

During optimization, you should see output like:

```
[I 2025-11-17 10:00:00,000] Trial 0 finished with value: 0.0234
[I 2025-11-17 10:00:01,000] Trial 1 finished with value: 0.0189
[I 2025-11-17 10:00:02,000] Trial 2 finished with value: inf  ← Bad trial (all filtered)
[I 2025-11-17 10:00:03,000] Trial 3 finished with value: 0.0156
```

Some `inf` values are normal (when bad parameter combinations filter everything). The optimizer learns from these and avoids similar combinations.

## Success Criteria

Optimization is working correctly if:
- ✓ Some trials finish with finite values (not all are inf)
- ✓ Best value decreases over time
- ✓ Filtered variants file is created
- ✓ Plots are generated

## Troubleshooting Checklist

- [x] **Data types correct**: IMF and other columns now converted to float
- [ ] **Column names match**: Config column names match data file columns
- [ ] **Required columns present**: SAMPLE, paternal_age, maternal_age, var_type, REF, ALT
- [ ] **Some variants pass filters**: Not all variants filtered out
- [ ] **Reference data valid**: Has required columns and variant types
- [ ] **Config file correct**: Valid YAML with correct optimization types
- [ ] **Sufficient trials**: At least 50-100 trials for meaningful optimization

## Next Steps

1. **Reinstall package**: `pip install -e . --force-reinstall --no-deps`
2. **Run optimization**: Use CLI command above with `-vv` flag
3. **Monitor progress**: Watch for finite trial values
4. **Check outputs**: Verify plots and filtered variants are created
5. **Review results**: Check filter_summary.txt for retention rates

If you continue to experience issues, share the output from running with `-vv` and I can help diagnose further.
