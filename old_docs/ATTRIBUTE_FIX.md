# Fix: 'PipelineConfig' object has no attribute 'optimization'

## Problem

After renaming the package, the code was still using the old attribute name `optimization` instead of `optimisation` in several places, causing the error:

```
Error during optimization: 'PipelineConfig' object has no attribute 'optimization'
```

## Root Cause

When we updated the configuration to use UK English spelling, we changed:
- YAML field: `optimization:` → `optimisation:`
- Class name: `OptimizationConfig` → `OptimisationConfig`
- Attribute: `config.optimization` → `config.optimisation`

However, the code in `pipeline.py` and `api.py` was still accessing `config.optimization`.

## Files Fixed

### 1. src/dnm_harmoniser/pipeline.py
Changed all references from `config.optimization` to `config.optimisation`:
- Line 125: `self.config.optimisation.variant_types`
- Line 139: `self.config.optimisation.regression_formula`
- Line 298: `self.config.optimisation.regression_formula`
- Line 303: `self.config.optimisation.regression_weights`
- Line 360: `self.config.optimisation.column_names`

### 2. src/dnm_harmoniser/api.py
Changed all references from `config.optimization` to `config.optimisation`:
- Line 308: `config.optimisation.variant_types`
- Line 309: `config.optimisation.regression_formula`

### 3. src/dnm_harmoniser/config.py
Already correct - uses `optimisation` throughout.

## Solution Applied

```bash
# Replace all .optimization with .optimisation in source files
find src/dnm_harmoniser -name "*.py" -type f -exec sed -i '' 's/\.optimization/.optimisation/g' {} +

# Reinstall package
pip install -e . --force-reinstall --no-deps
```

## Verification

After the fix:

```python
from dnm_harmoniser import PipelineConfig

config = PipelineConfig.from_yaml('/Users/kartikchundru/dnms/ukb/filter.yaml')

# All these now work correctly:
config.optimisation.variant_types
config.optimisation.regression_formula
config.optimisation.columns
config.optimisation.sample_id_column_data
```

✅ Configuration loads successfully
✅ All attributes accessible
✅ No AttributeError

## Summary

The attribute name throughout the codebase is now consistently `optimisation` (UK spelling) to match:
- The YAML field name
- The class name (`OptimisationConfig`)
- The attribute name on `PipelineConfig`

This ensures compatibility with configuration files using UK English spelling.
