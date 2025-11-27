"""Configuration models with validation using Pydantic.

Required Metadata Column Names:
    Specify the column names in your data for:
    - sample_id_column: Column name for sample IDs (default: "sample_id")
    - paternal_age_column: Column name for paternal age (default: "paternal_age")
    - maternal_age_column: Column name for maternal age (default: "maternal_age")
    - reference_column: Column name for reference allele (default: "ref")

Column Configuration:
    Each column can be configured with:
    - name: Column name in the dataframe
    - dtype: Data type (int, float, or str)
    - optimisation: Strategy to apply
        * "minimum": Set lower bound threshold
        * "maximum": Set upper bound threshold
        * "range": Keep values within specified range
        * "none": No thresholding (for other metadata columns)
    - range_constraint: Required only when optimisation="range"
        * min: Minimum allowed value
        * max: Maximum allowed value
    - linked_to: Optional name of another column to link with
        * Linked columns will always use the same threshold
        * Both columns must have the same optimisation type
        * Links must be bidirectional (both columns must reference each other)
    - variant_types: Optional list of variant types this column applies to
        * If None or omitted, applies to all variant types (SNV, Insertion, Deletion)
        * If specified, only applies to the listed variant types

Example:
    optimisation = OptimisationConfig(
        # Specify your column names
        sample_id_column="SampleID",
        paternal_age_column="Father_Age",
        maternal_age_column="Mother_Age",
        reference_column="REF",

        # Configure other columns for optimisation
        columns=[
            # Regular optimisation columns (apply to all variant types)
            ColumnConfig(name="quality_score", dtype="float", optimisation="minimum"),
            ColumnConfig(name="coverage", dtype="int", optimisation="maximum"),
            ColumnConfig(
                name="allele_balance",
                dtype="float",
                optimisation="range",
                range_constraint=RangeConstraint(min=0.25, max=0.75)
            ),

            # Variant-specific columns (only for SNVs)
            ColumnConfig(
                name="base_quality",
                dtype="float",
                optimisation="minimum",
                variant_types=["SNV"]
            ),

            # Indel-specific column (only for insertions and deletions)
            ColumnConfig(
                name="indel_length",
                dtype="int",
                optimisation="range",
                range_constraint=RangeConstraint(min=1, max=50),
                variant_types=["Insertion", "Deletion"]
            ),

            # Linked columns (share the same threshold)
            ColumnConfig(name="depth_father", dtype="int", optimisation="maximum", linked_to="depth_mother"),
            ColumnConfig(name="depth_mother", dtype="int", optimisation="maximum", linked_to="depth_father"),
        ]
    )
"""

from pathlib import Path
from typing import Dict, List, Literal, Optional, Any, Union
from pydantic import BaseModel, Field, field_validator, model_validator
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
    """Bayesian optimisation stage."""
    enabled: bool = True
    n_trials: int = Field(500, gt=0, description="Number of optimisation trials")
    sampler: Literal["tpe", "cmaes", "random"] = "tpe"
    pruner: Optional[Literal["successive_halving", "hyperband", "median"]] = "successive_halving"
    multivariate: bool = True


class RangeConstraint(BaseModel):
    """Range constraint for a column value."""
    min: Union[int, float]
    max: Union[int, float]

    @model_validator(mode='after')
    def validate_range(self):
        if self.min >= self.max:
            raise ValueError(f"min ({self.min}) must be less than max ({self.max})")
        return self


class SymmetricRangeConstraint(BaseModel):
    """
    Symmetric range constraint where upper = scale - lower.

    Used for allele balance, VAF, or other metrics that are symmetric around a midpoint.
    For example:
    - lower=0.25, scale=1.0 -> upper=0.75 (for 0-1 range)
    - lower=25, scale=100 -> upper=75 (for percentage range)
    """
    lower: Union[int, float] = Field(..., description="Lower bound (upper will be calculated as scale - lower)")
    scale: Union[int, float] = Field(1.0, description="Scale for symmetric constraint (1.0 for 0-1 range, 100 for percentage)")

    @model_validator(mode='after')
    def validate_symmetric_range(self):
        # Calculate upper bound
        upper = self.scale - self.lower

        # Validate that lower < midpoint < upper
        midpoint = self.scale / 2
        if self.lower >= midpoint:
            raise ValueError(
                f"lower ({self.lower}) must be less than midpoint ({midpoint}) for symmetric constraint"
            )

        # Validate bounds make sense
        if self.lower < 0:
            raise ValueError(f"lower ({self.lower}) must be non-negative")
        if upper > self.scale:
            raise ValueError(f"calculated upper ({upper}) exceeds scale ({self.scale})")

        return self

    @property
    def upper(self) -> Union[int, float]:
        """Calculate the upper bound from lower and scale."""
        return self.scale - self.lower

    @property
    def min(self) -> Union[int, float]:
        """Alias for lower bound (for compatibility with RangeConstraint)."""
        return self.lower

    @property
    def max(self) -> Union[int, float]:
        """Alias for upper bound (for compatibility with RangeConstraint)."""
        return self.upper


class ColumnConfig(BaseModel):
    """Configuration for a single column to optimise."""
    name: str = Field(..., description="Column name in the dataframe")
    dtype: Literal["int", "float", "str"] = Field(..., description="Data type of the column")
    optimisation: Literal["minimum", "maximum", "range", "none"] = Field(
        ...,
        description="Optimization strategy: minimum, maximum, range, or none (for metadata columns)"
    )
    range_constraint: Optional[Union[RangeConstraint, SymmetricRangeConstraint]] = Field(
        None,
        description="Range constraint (required if optimisation='range'). Use RangeConstraint for explicit min/max, or SymmetricRangeConstraint for symmetric bounds."
    )
    linked_to: Optional[str] = Field(
        None,
        description="Name of another column to link with (both will share the same threshold)"
    )
    variant_types: Optional[List[Literal["SNV", "Insertion", "Deletion"]]] = Field(
        None,
        description="Variant types this column applies to. If None, applies to all variant types."
    )

    @model_validator(mode='after')
    def validate_range_constraint(self):
        if self.optimisation == "range" and self.range_constraint is None:
            raise ValueError(f"range_constraint is required when optimisation='range' for column '{self.name}'")
        if self.optimisation != "range" and self.range_constraint is not None:
            raise ValueError(f"range_constraint should only be set when optimisation='range' for column '{self.name}'")
        if self.optimisation == "none" and self.linked_to is not None:
            raise ValueError(f"Metadata columns (optimisation='none') cannot be linked to other columns for column '{self.name}'")
        if self.variant_types is not None and len(self.variant_types) == 0:
            raise ValueError(f"variant_types must contain at least one variant type if specified for column '{self.name}'")
        return self


class OptimisationConfig(BaseModel):
    """Variant-specific optimisation settings."""
    variant_types: List[Literal["SNV", "Insertion", "Deletion"]] = ["SNV", "Insertion", "Deletion"]

    # Required metadata column names in input data
    sample_id_column_data: str = Field(
        "sample_id",
        description="Column name for sample ID in the input data"
    )
    paternal_age_column_data: str = Field(
        "paternal_age",
        description="Column name for paternal age in the input data"
    )
    maternal_age_column_data: str = Field(
        "maternal_age",
        description="Column name for maternal age in the input data"
    )
    reference_column_data: str = Field(
        "ref",
        description="Column name for reference allele in the input data"
    )
    alternate_column_data: str = Field(
        "alt",
        description="Column name for alternate allele in the input data"
    )

    # Required metadata column names in reference data
    sample_id_column_reference: str = Field(
        "sample_id",
        description="Column name for sample ID in the reference data"
    )
    paternal_age_column_reference: str = Field(
        "paternal_age",
        description="Column name for paternal age in the reference data"
    )
    maternal_age_column_reference: str = Field(
        "maternal_age",
        description="Column name for maternal age in the reference data"
    )
    reference_column_reference: str = Field(
        "ref",
        description="Column name for reference allele in the reference data"
    )
    alternate_column_reference: str = Field(
        "alt",
        description="Column name for alternate allele in the reference data"
    )

    columns: List[ColumnConfig] = Field(
        default_factory=list,
        description="Column configurations with types and optimisation criteria"
    )
    regression_formula: str = "dnm_count ~ paternal_age"
    regression_weights: List[float] = [1.0, 1.0]

    @model_validator(mode='after')
    def validate_linked_columns(self):
        """Validate that linked columns exist and form valid pairs."""
        column_names = {col.name for col in self.columns}
        linked_pairs = {}

        for col in self.columns:
            if col.linked_to is not None:
                # Check that the linked column exists
                if col.linked_to not in column_names:
                    raise ValueError(
                        f"Column '{col.name}' is linked to '{col.linked_to}', "
                        f"but '{col.linked_to}' is not defined in the configuration"
                    )

                # Check that linked columns have the same optimisation type
                linked_col = next(c for c in self.columns if c.name == col.linked_to)
                if col.optimisation != linked_col.optimisation:
                    raise ValueError(
                        f"Linked columns '{col.name}' and '{col.linked_to}' must have the same optimisation type. "
                        f"Got '{col.optimisation}' and '{linked_col.optimisation}'"
                    )

                # Check for bidirectional linking
                pair = tuple(sorted([col.name, col.linked_to]))
                if pair in linked_pairs:
                    if linked_pairs[pair] != col.optimisation:
                        raise ValueError(f"Inconsistent optimisation types for linked pair {pair}")
                else:
                    linked_pairs[pair] = col.optimisation

                # Verify bidirectional link
                if linked_col.linked_to != col.name:
                    raise ValueError(
                        f"Column '{col.name}' is linked to '{col.linked_to}', "
                        f"but '{col.linked_to}' is not linked back to '{col.name}'. "
                        f"Linked columns must reference each other."
                    )

        return self

    def get_required_metadata_columns_data(self) -> Dict[str, str]:
        """Get the required metadata column names for input data.

        Returns:
            Dictionary mapping field names to column names in the input data
        """
        return {
            "sample_id": self.sample_id_column_data,
            "paternal_age": self.paternal_age_column_data,
            "maternal_age": self.maternal_age_column_data,
            "reference": self.reference_column_data,
            "alternate": self.alternate_column_data,
        }

    def get_required_metadata_columns_reference(self) -> Dict[str, str]:
        """Get the required metadata column names for reference data.

        Returns:
            Dictionary mapping field names to column names in the reference data
        """
        return {
            "sample_id": self.sample_id_column_reference,
            "paternal_age": self.paternal_age_column_reference,
            "maternal_age": self.maternal_age_column_reference,
            "reference": self.reference_column_reference,
            "alternate": self.alternate_column_reference,
        }

    def get_metadata_columns(self) -> List[ColumnConfig]:
        """Get all additional metadata columns (optimisation='none')."""
        return [col for col in self.columns if col.optimisation == "none"]

    def get_optimisation_columns(self) -> List[ColumnConfig]:
        """Get all columns that will be optimised (excluding metadata)."""
        return [col for col in self.columns if col.optimisation != "none"]

    def get_columns_for_variant_type(self, variant_type: Literal["SNV", "Insertion", "Deletion"]) -> List[ColumnConfig]:
        """Get all columns that apply to a specific variant type.

        Args:
            variant_type: The variant type to filter by

        Returns:
            List of columns that apply to this variant type (either specifically or to all types)
        """
        return [
            col for col in self.columns
            if col.variant_types is None or variant_type in col.variant_types
        ]

    def get_optimisation_columns_for_variant_type(
        self,
        variant_type: Literal["SNV", "Insertion", "Deletion"]
    ) -> List[ColumnConfig]:
        """Get optimisation columns (excluding metadata) for a specific variant type.

        Args:
            variant_type: The variant type to filter by

        Returns:
            List of optimisation columns that apply to this variant type
        """
        return [
            col for col in self.get_optimisation_columns()
            if col.variant_types is None or variant_type in col.variant_types
        ]

    def get_linked_column_groups(self) -> List[List[ColumnConfig]]:
        """Get groups of linked columns."""
        seen = set()
        groups = []

        for col in self.columns:
            if col.linked_to is not None and col.name not in seen:
                # Find all columns in this linked group
                group = [col]
                seen.add(col.name)
                linked_col = next(c for c in self.columns if c.name == col.linked_to)
                group.append(linked_col)
                seen.add(linked_col.name)
                groups.append(group)

        return groups


class PipelineConfig(BaseModel):
    """Complete pipeline configuration."""
    stage1: Stage1Config = Field(default_factory=Stage1Config)
    stage2: Stage2Config = Field(default_factory=Stage2Config)
    stage3: Stage3Config = Field(default_factory=Stage3Config)
    optimisation: OptimisationConfig = Field(default_factory=OptimisationConfig)
    
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
        optimisation=OptimisationConfig(
            columns=[
                ColumnConfig(name="quality_score", dtype="float", optimisation="minimum"),
                ColumnConfig(name="coverage", dtype="int", optimisation="maximum"),
            ]
        ),
        max_workers=2
    ),
    "balanced": PipelineConfig(
        stage1=Stage1Config(n_trials=50),
        stage2=Stage2Config(min_dnm_count=5, max_dnm_count=300),
        stage3=Stage3Config(n_trials=500),
        optimisation=OptimisationConfig(
            columns=[
                ColumnConfig(name="quality_score", dtype="float", optimisation="minimum"),
                ColumnConfig(name="coverage", dtype="int", optimisation="maximum"),
                ColumnConfig(
                    name="allele_balance",
                    dtype="float",
                    optimisation="range",
                    range_constraint=RangeConstraint(min=0.25, max=0.75)
                ),
            ]
        ),
        max_workers=4
    ),
    "thorough": PipelineConfig(
        stage1=Stage1Config(n_trials=100),
        stage2=Stage2Config(min_dnm_count=10, max_dnm_count=200),
        stage3=Stage3Config(n_trials=1000, sampler="cmaes"),
        optimisation=OptimisationConfig(
            columns=[
                ColumnConfig(name="quality_score", dtype="float", optimisation="minimum"),
                ColumnConfig(name="coverage", dtype="int", optimisation="maximum"),
                ColumnConfig(
                    name="allele_balance",
                    dtype="float",
                    optimisation="range",
                    range_constraint=RangeConstraint(min=0.25, max=0.75)
                ),
            ]
        ),
        max_workers=8
    ),
    "noisy_data": PipelineConfig(
        stage1=Stage1Config(n_trials=50),
        stage2=Stage2Config(min_dnm_count=5, max_dnm_count=300),
        stage3=Stage3Config(n_trials=500, pruner="successive_halving"),
        optimisation=OptimisationConfig(
            columns=[
                ColumnConfig(name="quality_score", dtype="float", optimisation="minimum"),
                ColumnConfig(name="coverage", dtype="int", optimisation="maximum"),
            ]
        ),
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
