# Automatic Plotting Feature

## Overview

The optimization pipeline now automatically generates plots and saves filtered results after completing optimization. This eliminates the need to manually call plotting functions.

## What Gets Generated

When you run optimization with an `output_dir` specified, the pipeline automatically creates:

1. **Plots**: `optimization_results.png`
   - 2x3 grid layout
   - **Row 1**: Paternal age vs DNMs (SNV, Insertion, Deletion)
   - **Row 2**: Maternal age vs DNMs (SNV, Insertion, Deletion)
   - Shows both reference data (deCODE) and filtered input data (UKBB)

2. **Filtered Variants**: `filtered_variants.tsv`
   - All variants that pass the optimal filters
   - Includes a `variant_type` column for easy identification

3. **Summary Statistics**: `filter_summary.txt`
   - Original vs filtered variant counts
   - Retention percentages
   - Best parameters for each variant type

4. **Best Parameters**: `best_parameters.txt`
   - Human-readable format of optimal filtering parameters

5. **Parameters as YAML**: `optimal_params.yaml`
   - Machine-readable format for downstream analysis

6. **Optimization Summary**: `summary.txt`
   - Complete optimization results including scores

## Usage

### Simple API (Level 1)

```python
from dnm_harmoniser import optimize_filters

# Automatic plotting enabled by default when output_dir is specified
params = optimize_filters(
    data_path="ukbb_data.tsv",
    reference_path="decode_reference.tsv",
    output_dir="results/",  # Plots will be saved here
    preset="balanced",
    n_trials=100
)
```

### Advanced API (Level 2)

```python
from dnm_harmoniser import run_optimisation

result = run_optimisation(
    data="ukbb_data.tsv",
    reference="decode_reference.tsv",
    config_file="filter.yaml",
    output_dir="results/",      # Where to save plots and files
    generate_plots=True         # Enable/disable plotting (default: True)
)

# Access results
print(result.summary)
print(f"Best parameters: {result.best_params}")
```

### Disable Automatic Plotting

If you want to run optimization without generating plots:

```python
result = run_optimisation(
    data="ukbb_data.tsv",
    reference="decode_reference.tsv",
    config_file="filter.yaml",
    output_dir="results/",
    generate_plots=False  # Disable automatic plotting
)
```

### Manual Plotting

If you want to generate plots manually or customize them:

```python
from dnm_harmoniser import (
    run_optimisation,
    plot_optimization_results,
    save_parameters
)
from pathlib import Path

# Run optimization without automatic plotting
result = run_optimisation(
    data="ukbb_data.tsv",
    reference="decode_reference.tsv",
    config_file="filter.yaml",
    generate_plots=False
)

# Load your data
import pandas as pd
ukbb_df = pd.read_csv("ukbb_data.tsv", sep='\t')
decode_df = pd.read_csv("decode_reference.tsv", sep='\t')

# Generate plots manually with custom settings
plot_optimization_results(
    data_df=ukbb_df,
    reference_df=decode_df,
    best_params=result.best_params,
    output_dir=Path("custom_results/"),
    save_filtered=True
)

# Save parameters separately
save_parameters(
    result.best_params,
    output_dir=Path("custom_results/")
)
```

## Direct Pipeline Access (Level 3)

For full control over the pipeline:

```python
from dnm_harmoniser import PipelineConfig, OptimisationPipeline, VariantDataset
from pathlib import Path

# Load configuration
config = PipelineConfig.from_yaml("filter.yaml")

# Load data
data = VariantDataset.from_tsv(
    Path("ukbb_data.tsv"),
    sample_col=config.optimisation.sample_id_column_data,
    paternal_age_col=config.optimisation.paternal_age_column_data,
    maternal_age_col=config.optimisation.maternal_age_column_data,
    reference_col=config.optimisation.reference_column_data,
    alternate_col=config.optimisation.alternate_column_data
)

reference = VariantDataset.from_tsv(
    Path("decode_reference.tsv"),
    sample_col=config.optimisation.sample_id_column_reference,
    paternal_age_col=config.optimisation.paternal_age_column_reference,
    maternal_age_col=config.optimisation.maternal_age_column_reference,
    reference_col=config.optimisation.reference_column_reference,
    alternate_col=config.optimisation.alternate_column_reference
)

# Run pipeline with automatic plotting
pipeline = OptimisationPipeline(config)
result = pipeline.run(
    data,
    reference,
    output_dir=Path("results/"),
    generate_plots=True
)
```

## Output Directory Structure

After running optimization with `output_dir="results/"`:

```
results/
├── optimization_results.png     # 2x3 grid plot
├── filtered_variants.tsv        # Filtered variants
├── filter_summary.txt           # Summary statistics
├── best_parameters.txt          # Human-readable parameters
├── optimal_params.yaml          # Machine-readable parameters
└── summary.txt                  # Optimization summary
```

## Plot Details

### Layout
- **Figure size**: 24" × 14" (suitable for presentations)
- **Resolution**: 300 DPI (publication quality)
- **Grid**: 2 rows × 3 columns

### Data Shown
- **Reference**: deCODE data (blue regression line)
- **Filtered**: UKBB data after applying optimal filters (orange regression line)

### Axes
- **X-axis**: Parental age (years)
- **Y-axis**: Number of DNMs per person

### Variant Types
- **Column 1**: SNVs (Single Nucleotide Variants)
- **Column 2**: Insertions
- **Column 3**: Deletions

## Parameter Names

The filtered variants and plots use the new parameter naming scheme:

| Optimization Type | Parameter Format | Example |
|-------------------|------------------|---------|
| `minimum` | `max_{column}` | `max_DeNovoCNN_prob`, `max_DNM` |
| `maximum` | `min_{column}` | `min_child_coverage` |
| `range` | `min_{column}`, `max_{column}` | `min_VAF`, `max_VAF` |

See [PARAMETER_NAMING_GUIDE.md](PARAMETER_NAMING_GUIDE.md) for detailed explanation.

## Example: Complete Workflow

```python
from dnm_harmoniser import optimize_filters
from pathlib import Path

# Run optimization with automatic plotting
params = optimize_filters(
    data_path="/Users/kartikchundru/dnms/ukb/ukbb_data.tsv",
    reference_path="/Users/kartikchundru/dnms/decode/reference.tsv",
    output_dir="/Users/kartikchundru/dnms/results/",
    preset="balanced",
    n_trials=200
)

print("Optimization complete! Check results/ directory for:")
print("  - optimization_results.png")
print("  - filtered_variants.tsv")
print("  - filter_summary.txt")
print("  - best_parameters.txt")
```

## Key Features

1. **Automatic**: No need to manually call plotting functions
2. **Separate by Parent**: Paternal and maternal age plotted separately
3. **Separate by Variant Type**: SNV, Insertion, and Deletion plotted separately
4. **Complete Output**: Plots, filtered data, and summaries all saved automatically
5. **Optional**: Can be disabled with `generate_plots=False`
6. **Error Handling**: If plotting fails, optimization still completes successfully

## Notes

- Plotting only runs if `output_dir` is specified and `generate_plots=True`
- If plotting fails (e.g., missing columns), a warning is logged but optimization completes
- Filtered variants include all columns from the original data plus a `variant_type` column
- Summary statistics show retention percentages to help evaluate filter stringency
