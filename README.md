# DNM Harmoniser

Bayesian optimisation of _de novo_ variant (DNM) filtering thresholds.

The optimiser searches for per-variant-type (SNV, insertion, deletion) quality thresholds such that the
filtered callset's relationship between per-proband DNM count and paternal age matches a high-quality
reference (deCODE). It runs in three stages:

1. **Warmup**: a short exploratory search.
2. **Outlier removal**: drops probands whose filtered DNM count falls outside `[min_dnm_count, max_dnm_count]`.
3. **Full optimisation**: the main search on the cleaned sample set.

An optional **VAF quality penalty** adds a term to the loss that penalises filtered sets whose proband VAF
distribution departs from that expected for true heterozygous germline DNMs.

## Installation

```bash
pip install -e .

# Development tools (pytest, ruff, black, ...)
pip install -e ".[dev]"

# Optional memory usage logging
pip install -e ".[monitoring]"
```

Plotting libraries (matplotlib, seaborn) are core dependencies, so the default install produces plots
and filtered output.

## Quick start

```bash
# 1. (Optional) convert large TSVs to Parquet once, which avoids slow TSV reloads on every run
dnm-harmoniser convert data.tsv data.parquet
dnm-harmoniser convert reference.tsv reference.parquet

# 2. Check that the files load (uses default column names: SAMPLE, paternal_age, maternal_age, ref, alt)
dnm-harmoniser validate data.parquet --reference reference.parquet

# 3. Print the resolved configuration without running
dnm-harmoniser run data.parquet reference.parquet -c UKB_filter.yaml --dry-run

# 4. Run
dnm-harmoniser run data.parquet reference.parquet \
    -c UKB_filter.yaml \
    -o results/ \
    -w 32 -v
```

Data and reference files can be `.tsv`, `.tsv.gz` or `.parquet`; the format is detected from the extension.

### `run` options

| Option | Description |
|---|---|
| `-c, --config FILE` | YAML configuration (see below). Overrides `--preset`. |
| `-p, --preset NAME` | Built-in preset when no config is given (default `balanced`). |
| `-o, --output DIR` | Output directory (default `./results`). |
| `-n, --n-trials N` | Stage 3 trial count. |
| `--min-dnm / --max-dnm N` | Stage 2 outlier bounds. |
| `-w, --workers N` | Parallel workers. |
| `-s, --seed N` | Random seed (default 42). |
| `-O, --override KEY=VALUE` | Override any config field, repeatable, e.g. `-O stage3.n_trials=1000 -O optimisation.vaf_quality_weight=0.2`. |
| `-m, --memory-efficient` | Currently has no effect (both pipelines are the same class). |
| `--no-plots` | Skip plots **and** `filtered_variants.tsv`. |
| `--dry-run` | Show the configuration and exit. |
| `-v` / `-vv` | INFO / DEBUG logging. |

Other commands: `dnm-harmoniser show-presets`, `dnm-harmoniser init -o config.yaml` (writes a starting
config from a preset), `dnm-harmoniser --help`.

## Input data

Both the data and the reference need one row per variant with:

- proband / sample ID
- paternal age at birth
- maternal age at birth (not used by the default regression, so NA is fine)
- chromosome, position, reference allele, alternate allele

Column names can differ between the two files; map them in the config (`sample_id_column_data`,
`sample_id_column_reference`, etc.). Every quality metric you want to threshold must be listed under
`optimisation.columns` and exist in the data file.

Chromosome and position columns are detected automatically if they use a common name (`CHROM`, `chr`, `Chr`,
`Chromosome`, `POS`, `pos`, `Position`, ...). Otherwise, set `chromosome_column_data` / `position_column_data`
(and the `_reference` equivalents). This matters for Parquet input, which only loads the columns the config
needs: a coordinate column with an unrecognised name is left out of `filtered_variants.tsv`. Other columns not
named in the config are also left out with Parquet input, while TSV input keeps every column.

## Configuration

See `example_config.yaml` for a fully commented example, and `UKB_filter.yaml` / `AoU_filter.yaml` for
production configs.

### Filter columns

Each entry in `optimisation.columns` defines one threshold to optimise:

```yaml
columns:
  - name: MQ
    dtype: int                 # int | float | str
    optimisation: maximum      # keep values >= threshold
  - name: nparAADn0
    dtype: int
    optimisation: minimum      # keep values <= threshold
  - name: VAF
    dtype: float
    optimisation: range        # keep values within [lower, scale - lower]
    range_constraint:
      lower: 15
      scale: 100
  - name: Father coverage
    dtype: int
    optimisation: maximum
    linked_to: Mother coverage # share one threshold
  - name: IMF
    dtype: float
    optimisation: range
    range_constraint: {lower: 0.15, scale: 1}
    variant_types: [Insertion, Deletion]  # only applied to these types
```

- `range_constraint` with `min`/`max` bounds the search space; with `lower`/`scale` it defines a symmetric range.
- `computed: true` marks columns the loader derives from REF/ALT rather than reads from the input
  (`length`, `abs_length`); these need an explicit `range_constraint`.
- With `lower`/`scale`, the optimiser searches the lower threshold from `lower` up to `scale / 2` and sets
  the upper threshold to `scale - lower threshold`.
- `stage3.pruner` and `stage3.multivariate` are currently ignored (no pruning; TPE is always multivariate).

### Objective

The optimiser fits `dnm_count ~ paternal_age` per variant type and minimises the distance of the slope and
intercept from the reference fit. (`regression_formula` exists in the config but is currently ignored.)

```yaml
regression_weights: [1.0, 5.0]   # [intercept, slope] weights
loss_function: weighted_mse      # weighted_mse | absolute | huber | log_ratio | slope_priority | asymmetric | max_slope
```

Per-variant-type variants: `regression_weights_by_type`, `loss_function_by_type`.

### VAF quality penalty

Add to the `optimisation:` block:

```yaml
vaf_quality_metric: kl_divergence  # mean | proportion | skewness | kl_divergence; omit or null to disable
vaf_quality_weight: 0.1            # 0.1 = mild; 1.0 = equal weight to the regression loss
# vaf_quality_weight_by_type: {Insertion: 0.2, Deletion: 0.15}
vaf_column: VAF                    # default: proband_VAF
vaf_min: 0.25                      # KL divergence window
vaf_max: 0.75
```

> **The VAF column must be on a 0–1 scale.** The expected means, low-VAF thresholds and KL reference
> distribution are fixed on that scale (see `src/dnm_harmoniser/vaf_quality.py`). A percentage (0–100)
> column will give a meaningless penalty.

Metrics:

- `mean`: deviation of mean VAF from the expected value.
- `proportion`: excess fraction of low-VAF calls.
- `skewness`: excess left skew.
- `kl_divergence`: whole-shape divergence from a truncated normal.

All four use variant-type-specific expectations, because indels have natural alignment bias.

### Presets

Used when no `--config` is given (`dnm-harmoniser show-presets` for details):

| Preset | Warmup | Outlier removal | Main trials |
|---|---|---|---|
| `fast` | 30 | off | 100 |
| `balanced` | 50 | on | 500 |
| `thorough` | 100 | on | 1000 (CMA-ES) |
| `noisy_data` | 50 | on (5–300) | 500 |

## Outputs

Written to the `--output` directory:

| File | Contents |
|---|---|
| `optimal_params.yaml`, `best_parameters.txt` | Best thresholds per variant type |
| `filtered_variants.tsv` | Input variants passing the optimal filters (sample/age/allele/coordinate columns renamed to `SAMPLE`, `paternal_age`, `maternal_age`, `REF`, `ALT`, `CHROM`, `POS`) |
| `filter_summary.txt` | Counts before/after filtering per variant type |
| `summary.txt` | Optimisation summary |
| `config.yaml` | The exact configuration used |
| `optimization_results.png`, `regression_trajectory.png`, `optimization_summary.png`, `optuna_<type>/` | Regression, VAF and optimisation diagnostic plots |

## Python API

```python
from dnm_harmoniser import optimize_filters, run_optimisation, PipelineConfig

# Simple
params = optimize_filters("data.parquet", "reference.parquet", preset="balanced", n_trials=500)

# Full control
config = PipelineConfig.from_yaml("UKB_filter.yaml")
result = run_optimisation(data="data.parquet", reference="reference.parquet",
                          config=config, output_dir="results/")
print(result.summary)
```

## Helper scripts

Standalone pre/post-processing scripts in `src/dnm_harmoniser/`:

- `apply_chunks.py IN.tsv OUT.tsv REF_COL ALT_COL`: drops indels ≥20 bp and variants in homopolymers ≥8,
  reading in chunks (for very large inputs).
- `prox_filter.py IN.tsv OUT.tsv`: adds a `proximity_filter` column (0 = within 20 bp of a previously kept
  variant in the same sample). Expects columns `SAMPLE`, `Chromosome`, `Pos`, `ALT`.

## License

MIT License. See LICENSE for details.
