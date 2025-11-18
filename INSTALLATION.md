# Installation Guide

## Quick Start

1. **Install the package in development mode:**

```bash
pip install -e .
```

2. **Verify installation:**

```bash
python -c "import dnm_harmoniser; print('✓ Installation successful!')"
```

3. **Check CLI is available:**

```bash
dnm-harmoniser --help
```

## Package Structure

After installation, the package structure is:

```
src/
└── dnm_harmoniser/
    ├── __init__.py
    ├── api.py          # High-level API
    ├── cli.py          # Command-line interface
    ├── config.py       # Configuration models
    ├── data.py         # Data loading and validation
    └── pipeline.py     # Optimization pipeline
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

## Common Issues

### ModuleNotFoundError: No module named 'dnm_harmoniser'

**Solution:** Make sure you've run `pip install -e .` from the project root directory.

### Import errors after modifying code

**Solution:** The package is installed in editable mode, so changes to the code should be reflected immediately. If you're still having issues, try:

```bash
pip uninstall dnm-harmoniserr
pip install -e .
```

## Development Installation

For development with additional tools:

```bash
# Install with development dependencies
pip install -e ".[dev]"

# Install all optional dependencies
pip install -e ".[all]"
```

## Testing the Installation

Run this comprehensive test:

```python
from dnm_harmoniser.config import PipelineConfig
from pathlib import Path

# Load example configuration
config = PipelineConfig.from_yaml(Path('example_config.yaml'))

# Test all features
print("✓ Configuration loaded")
print(f"  - {len(config.optimization.columns)} columns configured")
print(f"  - {len(config.optimization.get_metadata_columns())} metadata columns")
print(f"  - {len(config.optimization.get_linked_column_groups())} linked groups")

# Test variant-specific columns
for vt in ['SNV', 'Insertion', 'Deletion']:
    cols = config.optimization.get_optimization_columns_for_variant_type(vt)
    print(f"  - {len(cols)} {vt} optimization columns")

print("✅ All tests passed!")
```
