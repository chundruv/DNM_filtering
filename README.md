# Variant Optimizer

Three-stage Bayesian optimization for genomic variant filtering with warmup, outlier removal, and full optimization.

## Features

- **Progressive API**: Simple function → Configuration object → Direct stage control
- **Three-stage optimization**: Warmup → Outlier removal → Full optimization
- **Performance**: Vectorized operations, caching, parallel processing
- **Reproducible**: Deterministic seeds, Docker support, pinned dependencies

## Installation

```bash
pip install -e .

# With visualization
pip install -e ".[viz]"

# For development
pip install -e ".[dev]"
```

## Quick Start

```python
from variant_optimizer import optimize_filters

# Simple usage
params = optimize_filters("data.tsv", "reference.tsv")

# With configuration
params = optimize_filters(
    "data.tsv", 
    "reference.tsv",
    preset="noisy_data",  # Optimized for noisy genomic data
    n_trials=500
)
```

## Command Line

```bash
# Basic optimization
variant-optimize data.tsv reference.tsv

# With outlier removal for noisy data
variant-optimize data.tsv reference.tsv \
    --preset noisy_data \
    --min-dnm 5 \
    --max-dnm 300

# Custom configuration
variant-optimize data.tsv reference.tsv \
    --config myconfig.yaml \
    --override stage3.n_trials=1000 \
    --workers 8
```

## Configuration Presets

- `fast`: Quick testing (30/0/100 trials)
- `balanced`: Default production (50/outlier/500 trials)
- `thorough`: Research-grade (100/outlier/1000 trials)
- `noisy_data`: For noisy genomics (50/outlier/500 trials + pruning)

## API Levels

### Level 1: Simple Function
```python
params = optimize_filters("data.tsv", "reference.tsv", preset="balanced")
```

### Level 2: Configuration Object
```python
from variant_optimizer import PipelineConfig, run_optimization

config = PipelineConfig(
    stage1={'n_trials': 100},
    stage2={'min_dnm_count': 10, 'max_dnm_count': 200},
    stage3={'n_trials': 1000, 'sampler': 'cmaes'}
)
result = run_optimization("data.tsv", "reference.tsv", config)
```

### Level 3: Direct Stage Control
```python
from variant_optimizer import StageRunner

runner = StageRunner()
warmup_params = runner.warmup(data, targets, n_trials=100)
clean_data = runner.remove_outliers(data, warmup_params, min_dnm=5, max_dnm=300)
final_params = runner.optimize(clean_data, targets, n_trials=1000)
```

## Performance Tips

- Use categorical dtypes for chromosomes, variant types
- Enable caching with `use_cache=True` (default)
- Use parallel workers: `--workers 8`
- Enable pruning for faster convergence

## Docker

```bash
docker build -t variant-optimizer .

docker run -v $(pwd)/data:/data \
           -v $(pwd)/results:/results \
           variant-optimizer \
           variant-optimize /data/variants.tsv /data/reference.tsv
```

## Citation

If you use this software, please cite:
```
@software{variant_optimizer,
  title = {Variant Optimizer: Three-stage Bayesian optimization for genomic variant filtering},
  author = {University of Exeter},
  year = {2024},
  url = {https://github.com/university-exeter/variant-optimizer}
}
```

## License

MIT License - See LICENSE file for details.
