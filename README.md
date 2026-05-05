# Variant Optimiser

Three-stage optimisation for _de novo_ variant filtering with warmup, outlier removal, and full optimisation.

## Installation

```bash
pip install -e .

# With visualisation
pip install -e ".[viz]"

# For development
pip install -e ".[dev]"
```

## Quick Start

```python
from dnm_harmoniser import optimize_filters

# Simple usage
params = optimize_filters("data.tsv", "reference.tsv")

# With configuration
params = optimize_filters(
    "data.tsv", 
    "reference.tsv",
    preset="balanced",
    n_trials=500
)
```

## Command Line

```bash
# Basic optimization
dnm-harmoniser data.tsv reference.tsv

# With outlier removal
dnm-harmoniser data.tsv reference.tsv \
    --preset balanced \
    --min-dnm 5 \
    --max-dnm 300

# Custom configuration
dnm-harmoniser data.tsv reference.tsv \
    --config myconfig.yaml \
    --override stage3.n_trials=1000 \
    --workers 8
```

## Configuration Presets

- `fast`: Quick testing (30/0/100 trials)
- `balanced`: Default production (50/outlier/500 trials)
- `thorough`: Research-grade (100/outlier/1000 trials)
- `noisy_data`: For noisy data (50/outlier/500 trials + pruning)

## Input data

_de novo_ variant data, and reference data (deCODE _de novo_ variant calls) are required in a .tsv format. Required columns are 
- Proband ID
- paternal age at birth
- maternal age at birth (although not used by default so NA is fine)
and variant
- chromosome
- position
- reference
- alternate 

All other genotype level quality metrics must be listed in the yaml config for optimisation.

See example yaml files for the config.

## License

MIT License - See LICENSE file for details.
