# Package Rename: variant-optimizer → dnm-harmoniser

## Summary

The package has been successfully renamed from `variant-optimizer` to `dnm-harmoniser`.

## Changes Made

### 1. Package Configuration
- **File**: `pyproject.toml`
  - Package name: `variant-optimizer` → `dnm-harmoniser`
  - CLI command: `variant-optimize` → `dnm-harmoniser`
  - Module name: `variant_optimizer` → `dnm_harmoniser`
  - URLs updated to reflect new repository name

### 2. Source Code
- **Directory renamed**: `src/variant_optimizer/` → `src/dnm_harmoniser/`
- **Module imports**: All working (using relative imports)
- **CLI name**: Updated in `cli.py` from `variant-optimize` to `dnm-harmoniser`
- **Docstrings**: Updated to reflect new package name

### 3. Documentation
All `.md` files updated:
- [INSTALLATION.md](INSTALLATION.md)
- [QUICK_REFERENCE.md](QUICK_REFERENCE.md)
- [COLUMN_CONFIG_GUIDE.md](COLUMN_CONFIG_GUIDE.md)
- [SYMMETRIC_RANGE_GUIDE.md](SYMMETRIC_RANGE_GUIDE.md)
- [SYMMETRIC_RANGE_UPDATE.md](SYMMETRIC_RANGE_UPDATE.md)
- [CHANGELOG.md](CHANGELOG.md)
- [LATEST_CHANGES.md](LATEST_CHANGES.md)
- [README.md](README.md)

### 4. Installation
- Old package uninstalled: `variant-optimizer`
- New package installed: `dnm-harmoniser`
- Installation mode: Editable (`pip install -e .`)

## New Usage

### Installation

```bash
pip install -e .
```

### Python API

```python
# Import the package
from dnm_harmoniser import PipelineConfig, ColumnConfig, SymmetricRangeConstraint

# Load configuration
config = PipelineConfig.from_yaml('example_config.yaml')

# Access configuration
print(config.optimisation.columns)
```

### CLI

```bash
# Show help
dnm-harmoniser --help

# Initialize configuration
dnm-harmoniser init --output config.yaml

# Show presets
dnm-harmoniser show-presets

# Validate data
dnm-harmoniser validate data.tsv --config config.yaml

# Run optimization
dnm-harmoniser run data.tsv reference.tsv --config config.yaml --output results/
```

## Migration Guide

If you have existing code using the old package name:

### Python Code

**Old:**
```python
from variant_optimizer import PipelineConfig
from variant_optimizer.config import ColumnConfig
```

**New:**
```python
from dnm_harmoniser import PipelineConfig
from dnm_harmoniser.config import ColumnConfig
```

### Command Line

**Old:**
```bash
variant-optimize run data.tsv --config config.yaml
```

**New:**
```bash
dnm-harmoniser run data.tsv reference.tsv --config config.yaml
```

### Configuration Files

No changes needed! YAML configuration files remain exactly the same.

## Verification

All functionality has been tested and verified:

✅ Package imports correctly
```bash
python -c "import dnm_harmoniser; print(dnm_harmoniser.__version__)"
# Output: 0.1.0
```

✅ CLI command works
```bash
dnm-harmoniser --help
```

✅ Configuration loading works
```python
from dnm_harmoniser import PipelineConfig
config = PipelineConfig.from_yaml('example_config.yaml')
```

✅ All features working:
- Symmetric range constraints
- Regular range constraints
- Linked columns
- Variant-specific columns
- UK English spelling (optimisation, minimum, maximum)

## Files Changed

### Core Package Files
- `pyproject.toml` - Package metadata and configuration
- `src/variant_optimizer/` → `src/dnm_harmoniser/` - Source directory renamed
- `src/dnm_harmoniser/__init__.py` - Package docstring updated
- `src/dnm_harmoniser/cli.py` - CLI name updated

### Documentation Files (9 files)
- All `.md` files updated with new package name
- CLI commands updated throughout documentation
- Import examples updated

### Configuration Files
- No changes required - YAML files work as-is

## Breaking Changes

⚠️ **For users of the old package:**

1. **Uninstall old package**:
   ```bash
   pip uninstall variant-optimizer
   ```

2. **Install new package**:
   ```bash
   pip install -e .
   ```

3. **Update imports** in your Python code:
   - `variant_optimizer` → `dnm_harmoniser`

4. **Update CLI commands**:
   - `variant-optimize` → `dnm-harmoniser`

## No Breaking Changes

✅ Configuration files (YAML) remain compatible
✅ All features work exactly the same
✅ API surface remains unchanged (just the package name changed)
