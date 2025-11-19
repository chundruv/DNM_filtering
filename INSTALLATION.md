# Installation Guide

## Quick Start

```bash
pip install -e .
```
## Usage

### Python API

```python
from dnm_harmoniser.config import PipelineConfig, ColumnConfig, RangeConstraint
from pathlib import Path

# Load configuration from YAML
config = PipelineConfig.from_yaml(Path('example_config.yaml'))

# Access configuration
print(f"Total columns: {len(config.optimization.columns)}")
print(f"Metadata columns: {len(config.optimization.get_metadata_columns())}")
print(f"Linked groups: {len(config.optimization.get_linked_column_groups())}")

# Get columns for specific variant type
snv_cols = config.optimization.get_optimization_columns_for_variant_type("SNV")
print(f"SNV optimization columns: {len(snv_cols)}")
```

### Command-Line Interface

```bash
# Show available commands
dnm-harmoniser --help

# Initialize a new configuration file
dnm-harmoniser init --output my_config.yaml

# Show available presets
dnm-harmoniser show-presets

# Validate your data files
dnm-harmoniser validate data.tsv --config my_config.yaml

# Run optimization
dnm-harmoniser run data.tsv --config my_config.yaml --output results/
```


