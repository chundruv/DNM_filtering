"""
DNM Harmoniser: Three-stage Bayesian optimisation for genomic variant filtering.

A fast, reproducible pipeline for optimizing variant filtering parameters
using warmup, outlier removal, and Bayesian optimisation.
"""

from .api import optimize_filters, run_optimisation
from .config import (
    PipelineConfig,
    Stage1Config,
    Stage2Config,
    Stage3Config,
    RangeConstraint,
    SymmetricRangeConstraint,
    ColumnConfig,
    OptimisationConfig
)
from .pipeline import OptimisationPipeline
from .data import VariantDataset
from .plotting import plot_optimization_results, apply_filters_from_params, save_parameters

__version__ = "0.1.0"

__all__ = [
    "optimize_filters",
    "run_optimisation",
    "PipelineConfig",
    "Stage1Config",
    "Stage2Config",
    "Stage3Config",
    "RangeConstraint",
    "SymmetricRangeConstraint",
    "ColumnConfig",
    "OptimisationConfig",
    "OptimisationPipeline",
    "VariantDataset",
    "plot_optimization_results",
    "apply_filters_from_params",
    "save_parameters"
]
