# Quick Start: Automatic Plotting

## One-Liner

```python
from dnm_harmoniser import optimize_filters

optimize_filters("data.tsv", "reference.tsv", output_dir="results/")
```

This automatically generates:
- `optimization_results.png` - 2×3 grid plot (paternal/maternal × SNV/Insertion/Deletion)
- `filtered_variants.tsv` - Filtered variants
- `filter_summary.txt` - Summary statistics
- `best_parameters.txt` - Optimal parameters

## Common Usage Patterns

### Basic Usage
```python
from dnm_harmoniser import optimize_filters

params = optimize_filters(
    data_path="ukbb_data.tsv",
    reference_path="decode_reference.tsv",
    output_dir="results/"  # Triggers automatic plotting
)
```

### With Configuration File
```python
from dnm_harmoniser import run_optimisation

result = run_optimisation(
    data="ukbb_data.tsv",
    reference="decode_reference.tsv",
    config_file="filter.yaml",
    output_dir="results/"
)
```

### Disable Automatic Plotting
```python
result = run_optimisation(
    data="ukbb_data.tsv",
    reference="decode_reference.tsv",
    config_file="filter.yaml",
    generate_plots=False  # No plots
)
```

### Manual Plotting
```python
from dnm_harmoniser import run_optimisation, plot_optimization_results
import pandas as pd

# Run without automatic plotting
result = run_optimisation(
    data="ukbb_data.tsv",
    reference="decode_reference.tsv",
    generate_plots=False
)

# Load data
ukbb_df = pd.read_csv("ukbb_data.tsv", sep='\t')
decode_df = pd.read_csv("decode_reference.tsv", sep='\t')

# Plot manually
plot_optimization_results(
    data_df=ukbb_df,
    reference_df=decode_df,
    best_params=result.best_params,
    output_dir="custom_results/"
)
```

## Output Files

```
results/
├── optimization_results.png     # 2×3 grid: (paternal, maternal) × (SNV, Ins, Del)
├── filtered_variants.tsv        # Variants passing optimal filters
├── filter_summary.txt           # Counts, percentages, parameters
├── best_parameters.txt          # Human-readable parameters
├── optimal_params.yaml          # Machine-readable parameters
└── summary.txt                  # Full optimization summary
```

## Plot Layout

```
           SNVs          Insertions     Deletions
         ┌─────────────┬─────────────┬─────────────┐
Paternal │             │             │             │
Age      │   Blue +    │   Blue +    │   Blue +    │
         │   Orange    │   Orange    │   Orange    │
         ├─────────────┼─────────────┼─────────────┤
Maternal │             │             │             │
Age      │   Blue +    │   Blue +    │   Blue +    │
         │   Orange    │   Orange    │   Orange    │
         └─────────────┴─────────────┴─────────────┘

Blue = Reference (deCODE)
Orange = Filtered Data (UKBB)
```

## Parameter Names Quick Reference

```python
# Your optimal parameters will look like:
{
    'SNV': {
        'max_DeNovoCNN_prob': 0.5,      # minimum opt → max_* param
        'max_DNM': 10.0,                 # minimum opt → max_* param
        'min_child_coverage': 10,        # maximum opt → min_* param
        'min_VAF': 30,                   # range opt → min_* and max_*
        'max_VAF': 70,
    },
    'Insertion': { ... },
    'Deletion': { ... }
}
```

**Rule of thumb:**
- **Minimize** column (lower is better) → use `max_*` param → filter OUT high values
- **Maximize** column (higher is better) → use `min_*` param → filter OUT low values
- **Range** column → use both `min_*` and `max_*` params

## Troubleshooting

| Problem | Solution |
|---------|----------|
| No plots generated | Add `output_dir="results/"` parameter |
| Plot is blank | Check that data has `var_type` column |
| Filtered file empty | Filters may be too stringent - check `filter_summary.txt` |
| Import error | Run `pip install -e .` to reinstall package |
| Missing columns | Update config file with correct column names |

## Key Configuration Columns

Your data must have these columns (or map to them in config):
- `SAMPLE` (or specify `sample_id_column_data` in config)
- `paternal_age` (or specify `paternal_age_column_data` in config)
- `maternal_age` (or specify `maternal_age_column_data` in config)
- `var_type` (values: "SNV", "Insertion", "Deletion")
- `REF`, `ALT` (reference and alternate alleles)

## Example Config Snippet

```yaml
optimisation:
  # Column names in your data files
  sample_id_column_data: SAMPLE      # Or "pid" if that's your column name
  paternal_age_column_data: paternal_age
  maternal_age_column_data: maternal_age
  reference_column_data: REF         # Or "Ref" if that's your column name
  alternate_column_data: ALT         # Or "Alt" if that's your column name

  columns:
    - name: DeNovoCNN_prob
      dtype: float
      optimisation: minimum  # → generates max_DeNovoCNN_prob parameter

    - name: child_coverage
      dtype: int
      optimisation: maximum  # → generates min_child_coverage parameter

    - name: VAF
      dtype: float
      optimisation: range    # → generates min_VAF and max_VAF parameters
      range_constraint:
        lower: 25
        scale: 100  # Symmetric: 25-75
```

## Documentation Files

- **Complete Guide**: [AUTOMATIC_PLOTTING_GUIDE.md](AUTOMATIC_PLOTTING_GUIDE.md)
- **Technical Details**: [PLOTTING_IMPLEMENTATION.md](PLOTTING_IMPLEMENTATION.md)
- **Implementation Status**: [IMPLEMENTATION_COMPLETE.md](IMPLEMENTATION_COMPLETE.md)
- **Parameter Naming**: [PARAMETER_NAMING_GUIDE.md](PARAMETER_NAMING_GUIDE.md)
- **Example Script**: [example_with_plotting.py](example_with_plotting.py)

## Checklist

Before running:
- [ ] Package installed: `pip install -e .`
- [ ] Data files exist and are readable
- [ ] Config file has correct column names
- [ ] Data has `var_type` column with values "SNV", "Insertion", "Deletion"
- [ ] Output directory is writable

After running:
- [ ] Check `optimization_results.png` for visual inspection
- [ ] Review `filter_summary.txt` for retention rates
- [ ] Verify `filtered_variants.tsv` has expected number of variants
- [ ] Check `best_parameters.txt` for parameter values

## Need Help?

1. Run the example: `python example_with_plotting.py`
2. Check the logs for warnings/errors
3. Verify your config matches the example above
4. Review the documentation files listed above
