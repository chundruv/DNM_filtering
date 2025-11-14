"""
Variant Optimizer: Three-stage Bayesian optimization for genomic variant filtering.

A fast, reproducible pipeline for optimizing variant filtering parameters
using warmup, outlier removal, and Bayesian optimization.
"""

from .api import optimize_filters, run_optimization
from .config import PipelineConfig, Stage1Config, Stage2Config, Stage3Config
from .pipeline import OptimizationPipeline
from .data import VariantDataset

__version__ = "0.1.0"

__all__ = [
    "optimize_filters",
    "run_optimization", 
    "PipelineConfig",
    "Stage1Config",
    "Stage2Config", 
    "Stage3Config",
    "OptimizationPipeline",
    "VariantDataset"
]
