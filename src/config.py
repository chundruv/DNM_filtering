"""Configuration models with validation using Pydantic."""

from pathlib import Path
from typing import Dict, List, Literal, Optional, Any
from pydantic import BaseModel, Field, field_validator
import yaml
import os


class Stage1Config(BaseModel):
    """Warmup stage configuration."""
    enabled: bool = True
    n_trials: int = Field(50, gt=0, description="Number of warmup trials")
    
    @field_validator('n_trials')
    @classmethod
    def validate_statistical_power(cls, v):
        if v < 30:
            raise ValueError("n_trials should be ≥30 for warmup statistical validity")
        return v


class Stage2Config(BaseModel):
    """Outlier removal stage."""
    enabled: bool = True
    min_dnm_count: Optional[int] = Field(5, gt=0, description="Minimum filtered DNMs")
    max_dnm_count: Optional[int] = Field(300, gt=0, description="Maximum filtered DNMs")
    
    @field_validator('max_dnm_count')
    @classmethod
    def validate_range(cls, v, info):
        if 'min_dnm_count' in info.data and info.data['min_dnm_count']:
            if v and v <= info.data['min_dnm_count']:
                raise ValueError("max_dnm_count must be > min_dnm_count")
        return v


class Stage3Config(BaseModel):
    """Bayesian optimization stage."""
    enabled: bool = True
    n_trials: int = Field(500, gt=0, description="Number of optimization trials")
    sampler: Literal["tpe", "cmaes", "random"] = "tpe"
    pruner: Optional[Literal["successive_halving", "hyperband", "median"]] = "successive_halving"
    multivariate: bool = True


class OptimizationConfig(BaseModel):
    """Variant-specific optimization settings."""
    variant_types: List[Literal["SNV", "Insertion", "Deletion"]] = ["SNV", "Insertion", "Deletion"]
    column_names: Optional[List[str]] = None  # None means use defaults
    regression_formula: str = "dnm_count ~ paternal_age + maternal_age"
    regression_weights: List[float] = [1.0, 1.0, 1.0]


class PipelineConfig(BaseModel):
    """Complete pipeline configuration."""
    stage1: Stage1Config = Field(default_factory=Stage1Config)
    stage2: Stage2Config = Field(default_factory=Stage2Config)
    stage3: Stage3Config = Field(default_factory=Stage3Config)
    optimization: OptimizationConfig = Field(default_factory=OptimizationConfig)
    
    # Performance settings
    max_workers: int = Field(4, gt=0, le=32)
    use_cache: bool = True
    cache_dir: Optional[Path] = None
    
    # Reproducibility
    seed: int = 42
    deterministic: bool = True
    
    @classmethod
    def from_yaml(cls, path: Path) -> "PipelineConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)
    
    def to_yaml(self, path: Path) -> None:
        """Save configuration to YAML file."""
        with open(path, 'w') as f:
            yaml.dump(self.model_dump(), f, default_flow_style=False)
    
    @classmethod
    def from_preset(cls, preset: str) -> "PipelineConfig":
        """Load a preset configuration."""
        return PRESETS[preset].copy()


# Preset configurations
PRESETS = {
    "fast": PipelineConfig(
        stage1=Stage1Config(n_trials=30),
        stage2=Stage2Config(enabled=False),
        stage3=Stage3Config(n_trials=100),
        max_workers=2
    ),
    "balanced": PipelineConfig(
        stage1=Stage1Config(n_trials=50),
        stage2=Stage2Config(min_dnm_count=5, max_dnm_count=300),
        stage3=Stage3Config(n_trials=500),
        max_workers=4
    ),
    "thorough": PipelineConfig(
        stage1=Stage1Config(n_trials=100),
        stage2=Stage2Config(min_dnm_count=10, max_dnm_count=200),
        stage3=Stage3Config(n_trials=1000, sampler="cmaes"),
        max_workers=8
    ),
    "noisy_data": PipelineConfig(
        stage1=Stage1Config(n_trials=50),
        stage2=Stage2Config(min_dnm_count=5, max_dnm_count=300),
        stage3=Stage3Config(n_trials=500, pruner="successive_halving"),
        max_workers=4,
        deterministic=True
    )
}


def load_config(
    config_file: Optional[Path] = None,
    preset: Optional[str] = None,
    cli_overrides: Optional[Dict[str, Any]] = None
) -> PipelineConfig:
    """
    Load configuration with proper precedence.
    
    Priority (highest to lowest):
    1. CLI overrides
    2. Config file
    3. Preset
    4. Defaults
    """
    if preset:
        config = PipelineConfig.from_preset(preset)
    else:
        config = PipelineConfig()
    
    if config_file and config_file.exists():
        file_config = PipelineConfig.from_yaml(config_file)
        config = config.model_copy(update=file_config.model_dump(exclude_unset=True))
    
    if cli_overrides:
        # Apply CLI overrides using dot notation
        config_dict = config.model_dump()
        for key, value in cli_overrides.items():
            keys = key.split('.')
            d = config_dict
            for k in keys[:-1]:
                d = d[k]
            d[keys[-1]] = value
        config = PipelineConfig(**config_dict)
    
    return config
