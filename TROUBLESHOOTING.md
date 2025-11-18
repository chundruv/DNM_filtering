# DNM Harmoniser Troubleshooting Guide

## Quick Diagnosis

If you're getting `inf` values for all trials, run the diagnostic tool:

```bash
python diagnose_inf_issue.py \
  /Users/kartikchundru/dnms/ukb/ukb_dnms_for_filtering.tsv \
  /Users/kartikchundru/dnms/decode_parages.txt \
  /Users/kartikchundru/dnms/ukb/filter.yaml
```

## Common Issues and Fixes

### Issue 1: All Trials Return `inf`

**Symptoms:**
```
[I 2025-11-17 15:25:58,631] Trial 0 finished with value: inf
[I 2025-11-17 15:25:59,123] Trial 1 finished with value: inf
[I 2025-11-17 15:26:00,456] Trial 2 finished with value: inf
```

**Causes:**
1. **Data type mismatch** - Columns stored as strings instead of numbers
2. **Column names mismatch** - Config column names don't match data file
3. **Too restrictive filtering** - Parameter combinations filter out all variants
4. **Missing required columns** - SAMPLE, paternal_age, maternal_age, var_type

**Fixes:**

1. **Check data types** (FIXED in latest version):
   - IMF, DNM, and other numeric columns now automatically converted
   - Diagnostic tool shows data types

2. **Verify column names**:
   ```python
   import pandas as pd
   df = pd.read_csv("your_data.tsv", sep='\t')
   print(df.columns.tolist())
   ```
   Make sure names match your `filter.yaml` config.

3. **Adjust parameter ranges**:
   - Run diagnostic to see suggested ranges
   - Modify `filter.yaml` if needed

### Issue 2: Linked Columns Have Different Values

**Symptoms:**
```python
{'min_father_coverage': 23, 'min_mother_coverage': 27}  # Should be same!
```

**Cause:** Bug in linked columns logic (FIXED in latest version)

**Fix:** Reinstall package:
```bash
pip install -e . --force-reinstall --no-deps
```

### Issue 3: Low Retention Rates

**Symptoms:**
```
Filtered: 62077 → 50 (0.08% retained)
```

**Causes:**
- Parameter ranges too restrictive
- Data quality issues
- Inappropriate optimization goals

**Fixes:**

1. **Review data distribution**:
   ```bash
   python diagnose_inf_issue.py ... | grep "DATA RANGES"
   ```

2. **Adjust optimization column ranges** in config:
   ```yaml
   columns:
     - name: child_coverage
       dtype: int
       optimisation: maximum
       # Add range_constraint if needed
       range_constraint:
         min: 10
         max: 100
   ```

3. **Start with fewer filters**:
   - Comment out some columns in config
   - Add them back incrementally

### Issue 4: No Plots Generated

**Symptoms:**
- Optimization completes but no `optimization_results.png`

**Causes:**
- `output_dir` not specified
- `generate_plots=False`
- Missing matplotlib/seaborn

**Fixes:**

1. **Specify output directory**:
   ```bash
   dnm-harmoniser run ... --output results/
   ```

2. **Check dependencies**:
   ```bash
   pip install matplotlib seaborn
   ```

3. **Enable plotting explicitly**:
   ```python
   result = run_optimisation(..., output_dir="results/", generate_plots=True)
   ```

### Issue 5: Memory Issues

**Symptoms:**
```
MemoryError: Unable to allocate ...
```

**Causes:**
- Large datasets
- Too many parallel workers
- Not enough RAM

**Fixes:**

1. **Reduce workers**:
   ```bash
   dnm-harmoniser run ... --workers 1
   ```

2. **Use sampling**:
   ```python
   # Sample data before optimization
   data_df = pd.read_csv("data.tsv", sep='\t')
   sampled = data_df.sample(n=10000, random_state=42)
   sampled.to_csv("sampled_data.tsv", sep='\t', index=False)
   ```

3. **Increase swap space** (OS-level)

## Diagnostic Checklist

Before running optimization:

- [ ] **Data file exists and is readable**
- [ ] **Reference file exists and is readable**
- [ ] **Config file is valid YAML**
- [ ] **All column names in config match data files**
- [ ] **Data has required columns**: SAMPLE, paternal_age, maternal_age, REF, ALT, var_type
- [ ] **Numeric columns are actually numeric** (not strings)
- [ ] **var_type column has values**: SNV, Insertion, Deletion
- [ ] **Package is installed**: `pip install -e .`

## Running Diagnostics

### 1. Quick Column Check

```bash
python -c "import pandas as pd; df = pd.read_csv('data.tsv', sep='\t'); print(df.columns.tolist()); print(df.dtypes)"
```

### 2. Full Diagnostic

```bash
python diagnose_inf_issue.py data.tsv reference.tsv config.yaml
```

Checks:
- Configuration validity
- Column availability
- Data types
- Value ranges
- Sample filtering
- Reference data

### 3. Test Run with Few Trials

```bash
dnm-harmoniser run data.tsv reference.tsv \
  --config filter.yaml \
  --n-trials 10 \
  --output test_results/ \
  -vv
```

Watch for:
- Finite values (not all inf)
- Decreasing best scores
- Reasonable retention rates

## Understanding Output

### Good Output:
```
[I 15:26:00] Trial 0 finished with value: 2.4356
[I 15:26:01] Trial 1 finished with value: 1.8923
[I 15:26:02] Trial 2 finished with value: inf       ← Some inf is OK
[I 15:26:03] Trial 3 finished with value: 1.7234
[I 15:26:04] Trial 4 finished with value: 1.6892  ← Improving!
```

### Bad Output:
```
[I 15:26:00] Trial 0 finished with value: inf
[I 15:26:01] Trial 1 finished with value: inf
[I 15:26:02] Trial 2 finished with value: inf
[I 15:26:03] Trial 3 finished with value: inf     ← All inf - something wrong!
```

## Expected Retention Rates

Based on your data:

| Variant Type | Typical Retention |
|--------------|-------------------|
| SNV | 30-70% |
| Insertion | 10-50% |
| Deletion | 10-50% |

If retention is:
- **< 1%**: Filters too stringent
- **> 95%**: Filters too lenient (not filtering enough)

## Parameter Interpretation

### Your Data Characteristics:

```
DeNovoCNN_prob: 0.901 - 1.000 (median: 0.990)
→ Very confident calls, most variants near 1.0

DNM: -0.335 - 0.000 (median: 0.000)
→ Most values exactly 0

MQ: 3 - 60 (median: 60)
→ Most have maximum mapping quality

nparAADn0: 0 - 2 (median: 0)
→ Most have 0 alternate reads in parents

child_coverage: 0 - 243 (median: 31)
→ Good spread, median coverage

VAF: 14 - 93 (median: 48)
→ Good spread around 50%
```

### Recommended Ranges:

```yaml
columns:
  - name: DeNovoCNN_prob
    dtype: float
    optimisation: minimum
    # Suggest 0.95-1.0 (most calls are confident)

  - name: child_coverage
    dtype: int
    optimisation: maximum
    # Suggest 20-50 (median is 31)

  - name: VAF
    dtype: float
    optimisation: range
    range_constraint:
      lower: 25    # Symmetric constraint
      scale: 100   # 25-75 range
```

## Getting Help

1. **Run diagnostic tool first**
2. **Check this troubleshooting guide**
3. **Review configuration examples**
4. **Check data format matches expectations**

## Files Reference

- **diagnose_inf_issue.py** - Main diagnostic tool
- **PARAMETER_NAMING_GUIDE.md** - Understand parameter naming
- **COLUMN_CONFIG_GUIDE.md** - Configure columns correctly
- **AUTOMATIC_PLOTTING_GUIDE.md** - Use automatic plotting
- **QUICK_START_PLOTTING.md** - Quick reference for plotting
- **README.md** - Main documentation

## Recent Fixes

### Version 0.1.0 (Current)

✅ **Fixed**: IMF and other columns converted to numeric automatically
✅ **Fixed**: Linked columns now work correctly (father_coverage = mother_coverage)
✅ **Fixed**: CLI now passes output_dir for automatic plotting
✅ **Added**: Comprehensive diagnostic tool
✅ **Added**: Automatic plotting integration

## Still Having Issues?

If none of the above helps:

1. Share the output of:
   ```bash
   python diagnose_inf_issue.py data.tsv reference.tsv config.yaml
   ```

2. Share the output of a test run with `-vv`:
   ```bash
   dnm-harmoniser run ... -vv > debug.log 2>&1
   ```

3. Check that your data format matches the expected format:
   - Tab-separated values
   - Headers in first row
   - No extra quotes or formatting
   - Numeric columns contain only numbers (or NaN)
