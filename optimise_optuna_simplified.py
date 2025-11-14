import pandas as pd
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf
import optuna
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, field

# --- Setup ---
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
sns.set_theme(style="whitegrid")


# --- Configuration Classes ---

@dataclass
class ColumnSpec:
    """Defines a single filtering column or group of linked columns."""
    column: str  # Can be a single column name or a descriptive name for a group
    param_type: str  # 'min', 'max', 'range', or 'weight'
    suggest_type: Optional[str] = None  # 'float' or 'int' - will be auto-detected if None
    bounds: Optional[Tuple[float, float]] = None  # Will be auto-detected if None
    transform: Optional[str] = None  # 'astype_float', etc.
    linked_columns: Optional[List[str]] = None  # If set, this threshold applies to all these columns
    
    def get_columns(self) -> List[str]:
        """Get list of columns this spec applies to."""
        if self.linked_columns:
            return self.linked_columns
        return [self.column]
    
    def get_param_name(self) -> str:
        """Generate parameter name for Optuna."""
        if self.param_type == 'min':
            return f'min_{self.column}'
        elif self.param_type == 'max':
            return f'max_{self.column}'
        elif self.param_type == 'range':
            return f'range_{self.column}'
        elif self.param_type == 'weight':
            return f'weight_{self.column}'
        return self.column


@dataclass
class FilterConfig:
    """Configuration for filtering and optimization."""
    variant_type: str
    columns: List[ColumnSpec]
    regression_formula: str = 'dnm_count ~ paternal_age + maternal_age'
    regression_weights: List[float] = field(default_factory=lambda: [1.0, 1.0, 1.0])
    n_trials: int = 500
    outlier_removal: Optional[Dict[str, Any]] = None  # Outlier removal configuration
    
    def remove_outliers(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Remove outliers from the data before optimization.
        
        Args:
            df: DataFrame to clean
            
        Returns:
            Cleaned DataFrame with outliers removed
        """
        if self.outlier_removal is None:
            return df
        
        df_clean = df.copy()
        initial_count = len(df_clean)
        
        method = self.outlier_removal.get('method', 'percentile')
        columns_to_clean = self.outlier_removal.get('columns', [])
        
        if not columns_to_clean:
            return df
        
        print(f"\n--- Removing outliers from {self.variant_type} data ---")
        print(f"Method: {method}")
        print(f"Initial variant count: {initial_count}")
        
        for col in columns_to_clean:
            if col not in df_clean.columns:
                print(f"Warning: Column '{col}' not found for outlier removal. Skipping.")
                continue
            
            col_data = pd.to_numeric(df_clean[col], errors='coerce')
            before_count = len(df_clean)
            
            if method == 'percentile':
                # Remove based on percentiles
                lower_percentile = self.outlier_removal.get('lower_percentile', 1)
                upper_percentile = self.outlier_removal.get('upper_percentile', 99)
                
                lower_bound = np.percentile(col_data.dropna(), lower_percentile)
                upper_bound = np.percentile(col_data.dropna(), upper_percentile)
                
                df_clean = df_clean[(col_data >= lower_bound) & (col_data <= upper_bound)]
                removed = before_count - len(df_clean)
                print(f"  {col}: Removed {removed} outliers "
                      f"(< {lower_bound:.2f} or > {upper_bound:.2f})")
                
            elif method == 'absolute':
                # Remove based on absolute values
                limits = self.outlier_removal.get('limits', {})
                if col in limits:
                    lower, upper = limits[col]
                    df_clean = df_clean[(col_data >= lower) & (col_data <= upper)]
                    removed = before_count - len(df_clean)
                    print(f"  {col}: Removed {removed} outliers "
                          f"(< {lower} or > {upper})")
                    
            elif method == 'iqr':
                # Remove based on IQR method
                q1 = col_data.quantile(0.25)
                q3 = col_data.quantile(0.75)
                iqr = q3 - q1
                multiplier = self.outlier_removal.get('iqr_multiplier', 1.5)
                
                lower_bound = q1 - multiplier * iqr
                upper_bound = q3 + multiplier * iqr
                
                df_clean = df_clean[(col_data >= lower_bound) & (col_data <= upper_bound)]
                removed = before_count - len(df_clean)
                print(f"  {col}: Removed {removed} outliers using IQR "
                      f"(< {lower_bound:.2f} or > {upper_bound:.2f})")
            
            elif method == 'zscore':
                # Remove based on z-score
                threshold = self.outlier_removal.get('zscore_threshold', 3)
                z_scores = np.abs((col_data - col_data.mean()) / col_data.std())
                df_clean = df_clean[z_scores <= threshold]
                removed = before_count - len(df_clean)
                print(f"  {col}: Removed {removed} outliers (|z-score| > {threshold})")
        
        final_count = len(df_clean)
        total_removed = initial_count - final_count
        pct_removed = (total_removed / initial_count) * 100
        
        print(f"Final variant count: {final_count}")
        print(f"Total removed: {total_removed} ({pct_removed:.1f}%)")
        
        return df_clean
    
    def auto_detect_bounds(self, df: pd.DataFrame) -> None:
        """Automatically detect bounds and types from data."""
        df_vartype = df[df['var_type'] == self.variant_type]
        
        for col_spec in self.columns:
            if col_spec.param_type == 'weight':
                # Weight parameters always have fixed bounds
                if col_spec.bounds is None:
                    col_spec.bounds = (0.0, 1.0)
                col_spec.suggest_type = 'float'
                continue
            
            # Get columns to analyze (single column or linked group)
            columns_to_check = col_spec.get_columns()
            
            # Collect data from all relevant columns
            all_col_data = []
            missing_cols = []
            
            for col in columns_to_check:
                if col not in df_vartype.columns:
                    missing_cols.append(col)
                    continue
                
                # Handle transforms
                col_data = df_vartype[col].copy()
                if col_spec.transform == 'astype_float':
                    col_data = col_data[col_data != '.']
                    col_data = pd.to_numeric(col_data, errors='coerce')
                else:
                    col_data = pd.to_numeric(col_data, errors='coerce')
                
                col_data = col_data.dropna()
                all_col_data.extend(col_data.tolist())
            
            if missing_cols:
                print(f"Warning: Column(s) {missing_cols} not found in data. Skipping.")
                if len(missing_cols) == len(columns_to_check):
                    continue
            
            if len(all_col_data) == 0:
                print(f"Warning: No valid data for column(s) {columns_to_check}. Skipping.")
                continue
            
            # Convert to numpy array for analysis
            combined_data = np.array(all_col_data)
            
            # Auto-detect suggest_type if not specified
            if col_spec.suggest_type is None:
                # Check if data is all integers
                if (combined_data % 1 == 0).all():
                    col_spec.suggest_type = 'int'
                else:
                    col_spec.suggest_type = 'float'
            
            # Auto-detect bounds if not specified
            if col_spec.bounds is None:
                data_min = combined_data.min()
                data_max = combined_data.max()
                
                if col_spec.param_type == 'min':
                    # For min filters, search from data_min to a bit above median
                    lower = data_min
                    upper = np.percentile(combined_data, 75)
                    if col_spec.suggest_type == 'int':
                        lower = int(np.floor(lower))
                        upper = int(np.ceil(upper))
                
                elif col_spec.param_type == 'max':
                    # For max filters, search from median to data_max
                    lower = np.percentile(combined_data, 25)
                    upper = data_max
                    if col_spec.suggest_type == 'int':
                        lower = int(np.floor(lower))
                        upper = int(np.ceil(upper))
                
                elif col_spec.param_type == 'range':
                    # For range filters (like VAF), determine symmetric range
                    if data_max <= 1.0:
                        # Data is 0-1, convert to percentage thinking
                        lower = 0.15
                        upper = 0.4
                    else:
                        # Data is likely 0-100
                        lower = 15
                        upper = 40
                
                col_spec.bounds = (lower, upper)
                
                # Print information about linked columns
                if col_spec.linked_columns:
                    print(f"Auto-detected bounds for {col_spec.column} (linked: {col_spec.linked_columns}) "
                          f"({col_spec.param_type}): {col_spec.bounds}, type: {col_spec.suggest_type}")
                else:
                    print(f"Auto-detected bounds for {col_spec.column} ({col_spec.param_type}): "
                          f"{col_spec.bounds}, type: {col_spec.suggest_type}")
    
    def get_optuna_params(self, trial) -> Dict[str, Any]:
        """Generate parameters for an Optuna trial."""
        params = {}
        for col_spec in self.columns:
            param_name = col_spec.get_param_name()
            if col_spec.suggest_type == 'float':
                params[param_name] = trial.suggest_float(param_name, col_spec.bounds[0], col_spec.bounds[1])
            elif col_spec.suggest_type == 'int':
                params[param_name] = trial.suggest_int(param_name, int(col_spec.bounds[0]), int(col_spec.bounds[1]))
        return params
    
    def apply_filters(self, df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
        """Apply filtering based on parameters."""
        filtered_df = df.copy()
        
        for col_spec in self.columns:
            if col_spec.param_type == 'weight':
                continue  # Weights don't filter data
            
            param_name = col_spec.get_param_name()
            param_value = params[param_name]
            
            # Get all columns this parameter applies to
            columns_to_filter = col_spec.get_columns()
            
            for col in columns_to_filter:
                if col not in filtered_df.columns:
                    continue
                
                # Apply any transforms
                if col_spec.transform == 'astype_float':
                    filtered_df = filtered_df[filtered_df[col] != '.'].copy()
                    filtered_df[col] = filtered_df[col].astype(float)
                
                # Apply filter
                if col_spec.param_type == 'min':
                    filtered_df = filtered_df[filtered_df[col] >= param_value]
                elif col_spec.param_type == 'max':
                    filtered_df = filtered_df[filtered_df[col] <= param_value]
                elif col_spec.param_type == 'range':
                    # For symmetric ranges (e.g., VAF between x and 100-x)
                    if filtered_df[col].max() <= 1.0:
                        # Data is 0-1
                        filtered_df = filtered_df[
                            (filtered_df[col] >= param_value) & 
                            (filtered_df[col] <= 1 - param_value)
                        ]
                    else:
                        # Data is 0-100
                        filtered_df = filtered_df[
                            (filtered_df[col] >= param_value) & 
                            (filtered_df[col] <= 100 - param_value)
                        ]
        
        return filtered_df
    
    def get_regression_weights(self, params: Dict[str, Any]) -> List[float]:
        """Get regression weights, including any optimized weights."""
        weights = self.regression_weights.copy()
        
        # Check if intercept weight is being optimized
        for i, col_spec in enumerate(self.columns):
            if col_spec.param_type == 'weight' and 'intercept' in col_spec.column.lower():
                param_name = col_spec.get_param_name()
                if param_name in params:
                    weights[0] = params[param_name]
                    break
        
        return weights


# --- Preset Configurations ---

def get_default_config(variant_type: str, column_names: Optional[List[str]] = None) -> FilterConfig:
    """
    Create a configuration with sensible defaults.
    
    Required columns (always included): paternal_age, maternal_age
    Default columns: paternal_age, maternal_age, VAF, child_DP, father_DP, mother_DP
    
    Note: father_DP and mother_DP are automatically linked to share the same threshold
    
    Args:
        variant_type: Type of variant (e.g., 'SNV', 'Insertion', 'Deletion')
        column_names: List of column names to filter on. If None, uses defaults.
                     Must include 'paternal_age' and 'maternal_age'.
    
    Returns:
        FilterConfig with automatic bound detection
    """
    # Required columns
    required = ['paternal_age', 'maternal_age']
    
    # Default columns if none specified
    if column_names is None:
        column_names = ['paternal_age', 'maternal_age', 'VAF', 'child_DP', 'father_DP', 'mother_DP']
    
    # Ensure required columns are present
    for req in required:
        if req not in column_names:
            raise ValueError(f"Required column '{req}' must be included in column_names")
    
    # Create column specifications
    columns = []
    processed_cols = set()
    
    for col in column_names:
        if col in processed_cols:
            continue
            
        if col in ['paternal_age', 'maternal_age']:
            # These are used in regression, not for filtering
            continue
        
        # Check if this is a parent depth column that should be linked
        if col in ['father_DP', 'mother_DP']:
            # Check if both parent depth columns are in the list
            if 'father_DP' in column_names and 'mother_DP' in column_names:
                # Create a single linked spec for both
                if 'father_DP' not in processed_cols:  # Only add once
                    columns.append(ColumnSpec(
                        'parent_DP',
                        'min',
                        linked_columns=['father_DP', 'mother_DP']
                    ))
                    processed_cols.add('father_DP')
                    processed_cols.add('mother_DP')
            else:
                # Only one parent DP specified, treat separately
                columns.append(ColumnSpec(col, 'min'))
                processed_cols.add(col)
        
        elif 'DP' in col or 'coverage' in col.lower() or 'depth' in col.lower():
            # Other depth/coverage columns: min filter
            columns.append(ColumnSpec(col, 'min'))
            processed_cols.add(col)
        
        elif 'VAF' in col or 'allele' in col.lower():
            # Allele frequency: range filter
            columns.append(ColumnSpec(col, 'range'))
            processed_cols.add(col)
        
        else:
            # Default: min filter
            columns.append(ColumnSpec(col, 'min'))
            processed_cols.add(col)
    
    # Add intercept weight optimization
    columns.append(ColumnSpec('intercept', 'weight'))
    
    return FilterConfig(
        variant_type=variant_type,
        columns=columns,
        regression_formula='dnm_count ~ paternal_age + maternal_age',
        regression_weights=[1.0, 1.0, 1.0],
        n_trials=500
    )


def get_custom_config(
    variant_type: str,
    column_specs: List[Dict[str, Any]],
    regression_formula: str = 'dnm_count ~ paternal_age + maternal_age',
    regression_weights: Optional[List[float]] = None,
    n_trials: int = 500
) -> FilterConfig:
    """
    Create a custom configuration with full control.
    
    Args:
        variant_type: Type of variant
        column_specs: List of dictionaries with column specifications:
            - column: str (required) - column name or descriptive name for linked group
            - param_type: str (required) - 'min', 'max', 'range', 'weight'
            - suggest_type: str (optional) - 'float' or 'int', auto-detected if not provided
            - bounds: tuple (optional) - (min, max), auto-detected if not provided
            - transform: str (optional) - 'astype_float', etc.
            - linked_columns: list (optional) - list of column names to apply same threshold to
        regression_formula: Formula for statsmodels regression
        regression_weights: Weights for regression parameters
        n_trials: Number of Optuna optimization trials
    
    Examples:
        # Single columns
        column_specs = [
            {'column': 'quality_score', 'param_type': 'min'},
            {'column': 'VAF', 'param_type': 'range', 'bounds': (20, 40)},
        ]
        
        # Linked columns (same threshold for multiple columns)
        column_specs = [
            {'column': 'parent_DP', 'param_type': 'min', 
             'linked_columns': ['father_DP', 'mother_DP']},
            {'column': 'child_DP', 'param_type': 'min'},
        ]
    """
    columns = [
        ColumnSpec(
            column=spec['column'],
            param_type=spec['param_type'],
            suggest_type=spec.get('suggest_type'),
            bounds=spec.get('bounds'),
            transform=spec.get('transform'),
            linked_columns=spec.get('linked_columns')
        )
        for spec in column_specs
    ]
    
    if regression_weights is None:
        regression_weights = [1.0] * (len(regression_formula.split('~')[1].split('+')) + 1)
    
    return FilterConfig(
        variant_type=variant_type,
        columns=columns,
        regression_formula=regression_formula,
        regression_weights=regression_weights,
        n_trials=n_trials
    )


# --- Data Loading and Preprocessing ---

def get_var_type(row):
    """Determines variant type based on REF and ALT allele lengths."""
    ref_len = len(str(row['REF']))
    alt_len = len(str(row['ALT']))
    if ref_len == 1 and alt_len == 1:
        return 'SNV'
    elif ref_len > alt_len:
        return 'Deletion'
    elif ref_len < alt_len:
        return 'Insertion'
    return 'Other'


def load_and_clean_data(filepath, sample_col='SAMPLE', required_cols=None, variant_type=None, max_length=20):
    """
    Loads and preprocesses the main data.
    
    Args:
        filepath: Path to data file
        sample_col: Name of sample ID column
        required_cols: List of required columns (besides SAMPLE, REF, ALT, paternal_age, maternal_age)
        variant_type: If provided and REF/ALT not present, creates var_type column with this value
        max_length: Maximum indel length
    """
    print("Loading data...")
    try:
        df = pd.read_csv(filepath, sep='\t', on_bad_lines='skip', low_memory=False)
        print("Data loaded successfully.")
    except Exception as e:
        print(f"Error loading file: {e}")
        return None

    # Clean quotes from column names and values
    df.columns = df.columns.str.replace('"', '')
    for col in df.select_dtypes(include=['object']).columns:
        if df[col].str.contains('"').any():
            df[col] = df[col].str.replace('"', '', regex=False)

    # Standardize sample column name
    if sample_col != 'SAMPLE' and sample_col in df.columns:
        df.rename(columns={sample_col: 'SAMPLE'}, inplace=True)
        sample_col = 'SAMPLE'

    # Convert common numeric columns
    numeric_cols = ['paternal_age', 'maternal_age', 'VAF']
    if required_cols:
        numeric_cols.extend([c for c in required_cols if c not in numeric_cols])
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Ensure REF and ALT are strings if they exist
    if 'Ref' in df.columns:
        df.rename(columns={'Ref': 'REF'}, inplace=True)
    elif 'Reference' in df.columns:
        df.rename(columns={'Reference': 'REF'}, inplace=True)
    if 'Alt' in df.columns:
        df.rename(columns={'Alt': 'ALT'}, inplace=True)
    elif 'Alternate' in df.columns:
        df.rename(columns={'Alternate': 'ALT'}, inplace=True)
    elif 'Variant' in df.columns:
        df.rename(columns={'Variant': 'ALT'}, inplace=True)

    # Drop rows with missing essential data
    essential_cols = ['SAMPLE', 'paternal_age', 'maternal_age']
    essential_cols = [c for c in essential_cols if c in df.columns]
    
    # Calculate midparage
    print("Calculating midparage and var_type...")
    df['midparage'] = (df['paternal_age'] + df['maternal_age']) / 2
    
    # Create variant type column
    if 'REF' in df.columns and 'ALT' in df.columns:
        # Calculate from REF/ALT
        df['var_type'] = df.apply(get_var_type, axis=1)
        df = df[df['var_type'] != 'Other']
    elif variant_type is not None:
        # Use provided variant type
        df['var_type'] = variant_type
        print(f"Created var_type column with value: {variant_type}")
    else:
        # Default to 'Unknown'
        df['var_type'] = 'Unknown'
        print("Warning: No REF/ALT columns and no variant_type specified. Using 'Unknown'.")
    
    df['length'] = df['ALT'].str.len() - df['REF'].str.len()

    df = df[np.abs(df['length']) < max_length]

    print(f"Loaded {len(df)} variants")
    return df


def load_reference_data(filepath, sample_col='SAMPLE', variant_type=None, max_length=20):
    """
    Loads and preprocesses the reference data.
    
    Args:
        filepath: Path to reference data file
        sample_col: Name of sample ID column in the file
        variant_type: If provided and REF/ALT not present, creates var_type column with this value
        max_length: Maximum indel length
    """
    print("\nLoading reference data...")
    try:
        df = pd.read_csv(filepath, sep='\t')
        
        # Standardize sample column name
        if 'pid' in df.columns:
            df.rename(columns={'pid': 'SAMPLE'}, inplace=True)
        elif sample_col != 'SAMPLE' and sample_col in df.columns:
            df.rename(columns={sample_col: 'SAMPLE'}, inplace=True)
        if 'Ref' in df.columns:
            df.rename(columns={'Ref': 'REF'}, inplace=True)
        elif 'Reference' in df.columns:
            df.rename(columns={'Reference': 'REF'}, inplace=True)
        if 'Alt' in df.columns:
            df.rename(columns={'Alt': 'ALT'}, inplace=True)
        elif 'Alternate' in df.columns:
            df.rename(columns={'Alternate': 'ALT'}, inplace=True)
        elif 'Variant' in df.columns:
            df.rename(columns={'Alternate': 'ALT'}, inplace=True)

        # Calculate midparage
        df['midparage'] = (df['paternal_age'] + df['maternal_age']) / 2
        
        # Create variant type column
        if 'REF' in df.columns and 'ALT' in df.columns:
            df['var_type'] = df.apply(get_var_type, axis=1)
            df = df[df['var_type'] != 'Other']
        elif variant_type is not None:
            df['var_type'] = variant_type
            print(f"Created var_type column with value: {variant_type}")
        else:
            df['var_type'] = 'Unknown'
            print("Warning: No REF/ALT columns and no variant_type specified. Using 'Unknown'.")
        
        df['length'] = df['ALT'].str.len() - df['REF'].str.len()

        df = df[np.abs(df['length']) < max_length]

        print("Reference data loaded successfully.")
        return df
    except Exception as e:
        print(f"Could not load reference file: {e}")
        return None


# --- Machine Learning Optimization ---

def remove_outlier_individuals(
    df: pd.DataFrame,
    configs: Dict[str, FilterConfig],
    warmup_params: Dict[str, Dict[str, Any]],
    min_dnm_count: Optional[int] = None,
    max_dnm_count: Optional[int] = None,
    sample_col: str = 'SAMPLE'
) -> pd.DataFrame:
    """
    Remove individuals based on TOTAL FILTERED DNM count across all variant types.
    
    This applies the warmup parameters to filter the data first, then counts DNMs,
    then removes outliers. This is more robust for noisy data.
    
    Args:
        df: DataFrame with variants (all variant types)
        configs: Dictionary of FilterConfig objects for each variant type
        warmup_params: Dictionary of warmup optimization parameters for each variant type
        min_dnm_count: Minimum total number of FILTERED DNMs (inclusive). If None, no lower bound.
        max_dnm_count: Maximum total number of FILTERED DNMs (inclusive). If None, no upper bound.
        sample_col: Name of sample ID column
    
    Returns:
        Filtered DataFrame with outlier individuals removed (across all variant types)
    """
    if min_dnm_count is None and max_dnm_count is None:
        return df  # No filtering needed
    
    print(f"\n{'='*60}")
    print("OUTLIER REMOVAL: Total Filtered DNMs Across All Variant Types")
    print(f"{'='*60}")
    print("\nApplying warmup filters and counting DNMs per individual...")
    
    # Apply warmup filters and count DNMs per individual for each variant type
    all_filtered_variants = []
    
    for var_type, config in configs.items():
        if var_type not in warmup_params:
            continue
        
        # Filter for this variant type
        df_vartype = df[df['var_type'] == var_type].copy()
        
        # Apply warmup filters
        filtered_df = config.apply_filters(df_vartype, warmup_params[var_type])
        
        print(f"  {var_type}: {len(df_vartype)} → {len(filtered_df)} variants after warmup filters")
        
        all_filtered_variants.append(filtered_df)
    
    # Combine all filtered variants
    if not all_filtered_variants:
        print("Warning: No filtered variants found!")
        return df
    
    df_filtered_all = pd.concat(all_filtered_variants, ignore_index=True)
    
    # Count TOTAL FILTERED DNMs per individual (across all variant types)
    total_dnm_counts = df_filtered_all.groupby(sample_col).size()
    
    print(f"\nTotal individuals before outlier removal: {len(total_dnm_counts)}")
    print(f"Filtered DNM count range: {total_dnm_counts.min():.0f} - {total_dnm_counts.max():.0f}")
    print(f"Mean filtered DNMs per individual: {total_dnm_counts.mean():.1f}")
    print(f"Median filtered DNMs per individual: {total_dnm_counts.median():.1f}")
    
    # Determine which samples to keep
    samples_to_keep = total_dnm_counts.index
    removed_low = 0
    removed_high = 0
    
    if min_dnm_count is not None:
        samples_low = total_dnm_counts[total_dnm_counts < min_dnm_count].index
        removed_low = len(samples_low)
        samples_to_keep = total_dnm_counts[total_dnm_counts >= min_dnm_count].index
        if removed_low > 0:
            print(f"\nRemoving {removed_low} individuals with < {min_dnm_count} filtered DNMs")
            counts_low = sorted(total_dnm_counts[samples_low].tolist())
            if len(counts_low) <= 20:
                print(f"  Their filtered DNM counts: {counts_low}")
            else:
                print(f"  Sample of their counts: {counts_low[:10]} ... {counts_low[-10:]}")
    
    if max_dnm_count is not None:
        if min_dnm_count is not None:
            samples_to_keep = total_dnm_counts[
                (total_dnm_counts >= min_dnm_count) & 
                (total_dnm_counts <= max_dnm_count)
            ].index
        else:
            samples_to_keep = total_dnm_counts[total_dnm_counts <= max_dnm_count].index
        
        samples_high = total_dnm_counts[total_dnm_counts > max_dnm_count].index
        removed_high = len(samples_high)
        if removed_high > 0:
            print(f"\nRemoving {removed_high} individuals with > {max_dnm_count} filtered DNMs")
            counts_high = sorted(total_dnm_counts[samples_high].tolist())
            if len(counts_high) <= 20:
                print(f"  Their filtered DNM counts: {counts_high}")
            else:
                print(f"  Sample of their counts: {counts_high[:10]} ... {counts_high[-10:]}")
    
    print(f"\n✓ Keeping {len(samples_to_keep)} individuals after outlier removal")
    if len(samples_to_keep) > 0:
        print(f"✓ Filtered DNM count range: {total_dnm_counts[samples_to_keep].min():.0f} - {total_dnm_counts[samples_to_keep].max():.0f}")
        print(f"✓ Filtered mean: {total_dnm_counts[samples_to_keep].mean():.1f}")
        print(f"✓ Filtered median: {total_dnm_counts[samples_to_keep].median():.1f}")
        
        # Show breakdown by variant type (after filtering)
        df_filtered = df[df[sample_col].isin(samples_to_keep)].copy()
        print(f"\nBreakdown by variant type (kept individuals, before full filters):")
        for vt in df_filtered['var_type'].unique():
            vt_total = len(df_filtered[df_filtered['var_type'] == vt])
            vt_per_ind = vt_total / len(samples_to_keep)
            print(f"  {vt}: {vt_total} variants ({vt_per_ind:.1f} per individual)")
    else:
        print("WARNING: No individuals remaining after outlier removal!")
    
    print(f"{'='*60}\n")
    
    # Filter original dataframe to keep only these samples
    return df[df[sample_col].isin(samples_to_keep)].copy()


def run_optimization_with_config(
    df: pd.DataFrame,
    config: FilterConfig,
    targets: np.ndarray,
    random_seed: int = 2025,
    sample_col: str = 'SAMPLE',
    remove_outliers: bool = True
) -> Optional[Dict[str, Any]]:
    """
    Runs Bayesian optimization using a configuration object.
    
    Args:
        df: DataFrame to optimize filtering on
        config: FilterConfig object defining the optimization
        targets: Target regression coefficients
        random_seed: Random seed for reproducibility
        sample_col: Name of sample ID column
        remove_outliers: Whether to remove outliers before optimization
    
    Returns:
        Dictionary of optimal parameters, or None if optimization fails
    """
    if df is None or df.empty:
        return None
    
    # Pre-filter for variant type
    df_vartype = df[df['var_type'] == config.variant_type].copy()
    
    # Remove outliers if requested
    if remove_outliers and config.outlier_removal is not None:
        df_vartype = config.remove_outliers(df_vartype)
        
        if df_vartype.empty:
            print("Warning: No data remaining after outlier removal!")
            return None
    
    parental_info = df_vartype[[sample_col, 'paternal_age', 'maternal_age']].drop_duplicates()
    
    def objective(trial):
        # Get parameters for this trial
        params = config.get_optuna_params(trial)
        
        # Apply filters
        try:
            filtered_df = config.apply_filters(df_vartype, params)
        except Exception as e:
            print(f"Error applying filters: {e}")
            return 1e6
        
        if filtered_df.empty:
            return 1e6

        # Count DNMs per sample
        dnm_counts = filtered_df.groupby(sample_col).size().rename('dnm_count').reset_index()
        regression_data = parental_info.merge(dnm_counts, on=sample_col, how='left').fillna(0)
        
        try:
            # Fit regression model
            model = smf.ols(config.regression_formula, data=regression_data).fit()
            model_params = model.params.values
            
            # Calculate weighted MSE
            squared_errors = (model_params - targets) ** 2
            weights = config.get_regression_weights(params)
            weighted_squared_errors = squared_errors * weights
            mse = np.mean(weighted_squared_errors)
        except Exception as e:
            print(f"Error in regression: {e}")
            return 1e6
        
        return mse

    print(f"\n--- Starting optimization for {config.variant_type} ---")
    sampler = optuna.samplers.TPESampler(seed=random_seed)
    study = optuna.create_study(direction="minimize", sampler=sampler)

    study.optimize(objective, n_trials=config.n_trials)

    print(f"\n--- Optimization for {config.variant_type} Finished ---")
    print(f"Best Score (MSE) found: {study.best_value:.6f}")
    print(f"Optimal filtering thresholds for {config.variant_type}:")
    for name, val in study.best_params.items():
        print(f"  - {name}: {val:.4f}" if isinstance(val, float) else f"  - {name}: {val}")
    
    return study.best_params


# --- Plotting ---

def plot_results(
    data_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    configs: Dict[str, FilterConfig],
    optimal_params: Dict[str, Dict[str, Any]],
    sample_col: str = 'SAMPLE'
):
    """Plot results for all variant types."""
    if data_df is None or reference_df is None or not optimal_params:
        print("\nSkipping plot due to missing data or optimization results.")
        return

    variant_types = list(optimal_params.keys())
    n_plots = len(variant_types)
    
    fig, axes = plt.subplots(1, n_plots, figsize=(8 * n_plots, 7), sharey=False)
    if n_plots == 1:
        axes = [axes]
    
    fig.suptitle('De Novo Mutations vs. Mid-Parental Age by Variant Type', 
                 fontsize=20, fontweight='bold')

    for i, var_type in enumerate(variant_types):
        ax = axes[i]
        config = configs[var_type]
        best_params = optimal_params[var_type]
        
        # Filter data with optimal parameters
        data_vartype = data_df[data_df['var_type'] == var_type]
        filtered_data_df = config.apply_filters(data_vartype, best_params)
        
        # Prepare reference data
        reference_subset = reference_df[reference_df['var_type'] == var_type]
        reference_counts = reference_subset.groupby(sample_col).size().rename('dnm_count')
        reference_ages = reference_subset[[sample_col, 'midparage']].drop_duplicates().set_index(sample_col)
        plot_data_reference = reference_ages.join(reference_counts, how='left').fillna(0)
        
        # Prepare data
        data_counts = filtered_data_df.groupby(sample_col).size().rename('dnm_count')
        data_ages = filtered_data_df[[sample_col, 'midparage']].drop_duplicates().set_index(sample_col)
        plot_data_data = data_ages.join(data_counts, how='left').fillna(0)

        # Plotting
        sns.regplot(data=plot_data_reference, x='midparage', y='dnm_count',
                    label='Reference', color='royalblue', ax=ax,
                    scatter_kws={'alpha': 0.6, 's': 40, 'edgecolor': 'w'},
                    line_kws={'linewidth': 2})
        sns.regplot(data=plot_data_data, x='midparage', y='dnm_count',
                    label='Data (Filtered)', color='darkorange', ax=ax,
                    scatter_kws={'alpha': 0.7, 's': 40, 'edgecolor': 'w'},
                    line_kws={'linewidth': 2})
        
        ax.set_title(f'{var_type}', fontsize=16)
        ax.set_xlabel('Mid-Parental Age', fontsize=14)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend()

    axes[0].set_ylabel('Number of DNMs per Person', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('dnm_out/dnm_optimization_results.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved to dnm_out/dnm_optimization_results.png")
    plt.show()


# --- Main Execution ---


def main(
    data_path: str,
    reference_path: str,
    configs: Optional[Dict[str, FilterConfig]] = None,
    column_names: Optional[List[str]] = None,
    variant_types: Optional[List[str]] = None,
    random_seed: int = 2025,
    sample_col_data: str = 'SAMPLE',
    sample_col_reference: str = 'SAMPLE',
    min_dnm_count: Optional[int] = None,
    max_dnm_count: Optional[int] = None,
    warmup_trials: Optional[int] = 100,
    refine_after_outlier_removal: bool = True
):
    """
    Main optimization pipeline with warmup and outlier removal support.
    
    THREE-STAGE PROCESS FOR NOISY DATA:
    1. Warmup: Quick optimization (warmup_trials) to get initial filters
    2. Outlier Removal: Apply warmup filters, count DNMs, remove outliers
    3. Full Optimization: Run full optimization (n_trials) on cleaned data
    
    Args:
        data_path: Path to main data file
        reference_path: Path to reference data file
        configs: Dictionary mapping variant types to FilterConfig objects.
                 If None, uses automatic configuration.
        column_names: List of column names to use for filtering (if configs is None).
                     Default: ['paternal_age', 'maternal_age', 'VAF', 'child_DP', 'father_DP', 'mother_DP']
                     Must include 'paternal_age' and 'maternal_age'.
        variant_types: List of variant types to process. Default: ['SNV', 'Insertion', 'Deletion']
        random_seed: Random seed for reproducibility
        sample_col_data: Name of sample ID column in data file
        sample_col_reference: Name of sample ID column in reference file
        min_dnm_count: Minimum TOTAL number of FILTERED DNMs per individual (inclusive).
                      After applying warmup filters and counting DNMs, individuals with fewer
                      total filtered DNMs are removed.
        max_dnm_count: Maximum TOTAL number of FILTERED DNMs per individual (inclusive).
                      After applying warmup filters and counting DNMs, individuals with more
                      total filtered DNMs are removed.
        warmup_trials: Number of trials for warmup optimization. Default: 100.
                      Set to None or 0 to skip warmup (outlier removal happens on raw counts).
                      Warmup optimization is faster and provides initial filters for outlier detection.
        refine_after_outlier_removal: DEPRECATED - kept for backwards compatibility.
                                     When warmup_trials is set, full optimization always runs after
                                     outlier removal.
    
    Returns:
        Dictionary of optimal parameters for each variant type.
    
    Examples:
        # Basic usage - no outlier removal
        results = main('data.tsv', 'reference.tsv')
        
        # With warmup and outlier removal (RECOMMENDED for noisy data)
        results = main(
            'data.tsv', 'reference.tsv',
            min_dnm_count=5,      # Remove individuals with <5 filtered DNMs
            max_dnm_count=300,    # Remove individuals with >300 filtered DNMs
            warmup_trials=100     # Quick warmup to get initial filters
        )
        # Process:
        # 1. Warmup: 100 trials per variant type
        # 2. Apply warmup filters and count DNMs
        # 3. Remove individuals with <5 or >300 filtered DNMs
        # 4. Full optimization: 500 trials (from config.n_trials) on cleaned data
        
        # Without warmup (outlier removal on raw counts)
        results = main(
            'data.tsv', 'reference.tsv',
            min_dnm_count=10,
            max_dnm_count=500,
            warmup_trials=None    # No warmup, use raw DNM counts
        )
    """
    np.random.seed(random_seed)
    
    # Load data
    data_df_all = load_and_clean_data(data_path, sample_col=sample_col_data, required_cols=column_names)
    data_df = data_df_all.copy()
    data_df.dropna(subset=['paternal_age', 'maternal_age'], inplace=True)

    reference_df = load_reference_data(reference_path, sample_col=sample_col_reference)
    
    if data_df is None or reference_df is None:
        print("Error: Could not load required data files.")
        return None
    
    # Standardize sample column to 'SAMPLE' for internal use
    sample_col = 'SAMPLE'
    
    # Determine variant types to process
    if variant_types is None:
        variant_types = ['SNV', 'Insertion', 'Deletion']
    
    # Create configs if not provided
    if configs is None:
        print("\nCreating automatic configurations...")
        configs = {}
        for vt in variant_types:
            configs[vt] = get_default_config(vt, column_names)
            print(f"Auto-detecting bounds for {vt}...")
            configs[vt].auto_detect_bounds(data_df)
    
    # Determine if we're doing three-stage optimization (warmup + outlier removal + full)
    doing_warmup = warmup_trials is not None and warmup_trials > 0
    doing_outlier_removal = (min_dnm_count is not None or max_dnm_count is not None)
    
    if doing_warmup and doing_outlier_removal:
        print(f"\n{'='*60}")
        print("THREE-STAGE OPTIMIZATION FOR NOISY DATA")
        print(f"{'='*60}")
        print(f"Stage 1: Warmup optimization ({warmup_trials} trials per variant)")
        print(f"Stage 2: Outlier removal (< {min_dnm_count} or > {max_dnm_count} filtered DNMs)")
        print(f"Stage 3: Full optimization on cleaned data")
        print(f"{'='*60}\n")
    
    data_df_original = data_df.copy()
    all_optimal_params = {}
    all_warmup_params = {}
    
    # STAGE 1: WARMUP OPTIMIZATION (if requested)
    if doing_warmup and doing_outlier_removal:
        print(f"\n{'#'*60}")
        print("STAGE 1: WARMUP OPTIMIZATION")
        print(f"{'#'*60}\n")
        
        for var_type in variant_types:
            if var_type not in configs:
                print(f"No configuration for {var_type}. Skipping.")
                continue
            
            config = configs[var_type]
            
            print(f"\n{'='*60}")
            print(f"Warmup Optimization: {var_type}")
            print(f"{'='*60}")
            
            # Calculate targets from reference data
            reference_subset = reference_df[reference_df['var_type'] == var_type]
            reference_counts = reference_subset.groupby(sample_col).size().rename('dnm_count')
            reference_ages = reference_subset[[sample_col, 'paternal_age', 'maternal_age']].drop_duplicates().set_index(sample_col)
            reference_regression_data = reference_ages.join(reference_counts, how='left').fillna(0)
            
            if reference_regression_data.shape[0] <= 2:
                print(f"Not enough reference data for {var_type}. Skipping.")
                continue
            
            try:
                model_reference = smf.ols(config.regression_formula, data=reference_regression_data).fit()
                targets = model_reference.params.values
                
                # Create warmup config with fewer trials
                warmup_config = FilterConfig(
                    variant_type=config.variant_type,
                    columns=config.columns,
                    regression_formula=config.regression_formula,
                    regression_weights=config.regression_weights,
                    n_trials=warmup_trials
                )
                
                # Run warmup optimization
                warmup_params = run_optimization_with_config(
                    data_df, warmup_config, targets, random_seed, sample_col
                )
                
                if warmup_params:
                    all_warmup_params[var_type] = warmup_params
                    print(f"✓ Warmup complete for {var_type}")
                else:
                    print(f"✗ Warmup failed for {var_type}")
                
            except Exception as e:
                print(f"Error in warmup for {var_type}: {e}")
                continue
        
        # STAGE 2: OUTLIER REMOVAL
        if all_warmup_params:
            print(f"\n{'#'*60}")
            print("STAGE 2: OUTLIER REMOVAL")
            print(f"{'#'*60}")
            
            data_df = remove_outlier_individuals(
                data_df,
                configs,
                all_warmup_params,
                min_dnm_count,
                max_dnm_count,
                sample_col
            )
            
            print(f"✓ Outlier removal complete")
            print(f"✓ Proceeding with {len(data_df[sample_col].unique())} individuals\n")
        else:
            print("\nWARNING: No warmup parameters available. Skipping outlier removal.")
            print("Proceeding with full dataset.\n")
        
        # STAGE 3: FULL OPTIMIZATION
        print(f"\n{'#'*60}")
        print("STAGE 3: FULL OPTIMIZATION ON CLEANED DATA")
        print(f"{'#'*60}\n")
    
    elif doing_outlier_removal and not doing_warmup:
        # Outlier removal without warmup (use raw counts)
        print(f"\n{'='*60}")
        print("OUTLIER REMOVAL (No Warmup)")
        print(f"{'='*60}\n")
        print("WARNING: Removing outliers based on RAW DNM counts (no filters applied).")
        print("Consider using warmup_trials parameter for better outlier detection.\n")
        
        total_dnm_counts = data_df.groupby(sample_col).size()
        samples_to_keep = total_dnm_counts.index
        
        if min_dnm_count is not None:
            samples_to_keep = total_dnm_counts[total_dnm_counts >= min_dnm_count].index
        if max_dnm_count is not None:
            if min_dnm_count is not None:
                samples_to_keep = total_dnm_counts[
                    (total_dnm_counts >= min_dnm_count) & 
                    (total_dnm_counts <= max_dnm_count)
                ].index
            else:
                samples_to_keep = total_dnm_counts[total_dnm_counts <= max_dnm_count].index
        
        removed = len(total_dnm_counts) - len(samples_to_keep)
        print(f"Removed {removed} individuals")
        print(f"Keeping {len(samples_to_keep)} individuals\n")
        
        data_df = data_df[data_df[sample_col].isin(samples_to_keep)].copy()
    
    # Run final optimization for each variant type
    for var_type in variant_types:
        if var_type not in configs:
            print(f"No configuration for {var_type}. Skipping.")
            continue
        
        config = configs[var_type]
        
        print(f"\n{'='*60}")
        print(f"Final Optimization: {var_type}")
        print(f"{'='*60}")
        
        # Calculate targets from reference data
        reference_subset = reference_df[reference_df['var_type'] == var_type]
        reference_counts = reference_subset.groupby(sample_col).size().rename('dnm_count')
        reference_ages = reference_subset[[sample_col, 'paternal_age', 'maternal_age']].drop_duplicates().set_index(sample_col)
        reference_regression_data = reference_ages.join(reference_counts, how='left').fillna(0)
        
        if reference_regression_data.shape[0] <= 2:
            print(f"Not enough reference data for {var_type}. Skipping.")
            continue
        
        try:
            model_reference = smf.ols(config.regression_formula, data=reference_regression_data).fit()
            targets = model_reference.params.values
            print(f"\nTarget coefficients for {var_type}:")
            for name, val in zip(model_reference.params.index, targets):
                print(f"  - {name}: {val:.4f}")
            
            # Run full optimization
            optimal_params = run_optimization_with_config(
                data_df, config, targets, random_seed, sample_col
            )
            
            if optimal_params:
                all_optimal_params[var_type] = optimal_params
                
                # Show comparison with warmup if available
                if var_type in all_warmup_params:
                    print(f"\n--- Comparison: Warmup vs Final Parameters ---")
                    for param_name in optimal_params.keys():
                        warmup_val = all_warmup_params[var_type].get(param_name, 'N/A')
                        final_val = optimal_params[param_name]
                        if isinstance(final_val, float) and isinstance(warmup_val, float):
                            print(f"  {param_name}: {warmup_val:.4f} → {final_val:.4f}")
                        else:
                            print(f"  {param_name}: {warmup_val} → {final_val}")
            else:
                print(f"Final optimization failed for {var_type}. Skipping.")
            
        except Exception as e:
            print(f"Error processing {var_type}: {e}")
            import traceback
            traceback.print_exc()
            continue

    os.mkdir('dnm_out/')
    # Plot results
    if all_optimal_params:
        plot_results(data_df, reference_df, configs, all_optimal_params, sample_col)
    
    config.apply_filters(data_df_all[data_df_all['var_type']=='SNV'], all_optimal_params['SNV']).to_csv('dnm_out/final_snvs.csv', index=False)
    config.apply_filters(data_df_all[data_df_all['var_type']=='Insertion'], all_optimal_params['SNV']).to_csv('dnm_out/final_insertion.csv', index=False)
    config.apply_filters(data_df_all[data_df_all['var_type']=='Deletion'], all_optimal_params['SNV']).to_csv('dnm_out/final_deletions.csv', index=False)

    return all_optimal_params


if __name__ == "__main__":
    RANDOM_SEED = 2025
    
    # Example 1: Use defaults (simplest)
    # Note: father_DP and mother_DP automatically linked by default!
    optimal_params = main(
        data_path='data_snv_dnms_for_filtering.tsv',
        reference_path='reference_parages.txt',
        random_seed=RANDOM_SEED
    )
    
    # Example 2: With outlier removal (removes based on TOTAL DNMs across all variant types)
    # optimal_params = main(
    #     data_path='data_snv_dnms_for_filtering.tsv',
    #     reference_path='reference_parages.txt',
    #     min_dnm_count=5,    # Remove individuals with < 5 TOTAL DNMs
    #     max_dnm_count=300,  # Remove individuals with > 300 TOTAL DNMs
    #     random_seed=RANDOM_SEED
    # )
    
    # Example 3: Custom sample columns and outlier removal
    # optimal_params = main(
    #     data_path='data_snv_dnms_for_filtering.tsv',
    #     reference_path='reference_parages.txt',
    #     sample_col_data='patient_id',      # Sample column in data file
    #     sample_col_reference='subject_id',  # Sample column in reference file
    #     min_dnm_count=10,
    #     max_dnm_count=250,
    #     random_seed=RANDOM_SEED
    # )
    
    # Example 4: Manual configuration with linked columns and outlier removal
    # custom_configs = {
    #     'SNV': get_custom_config(
    #         variant_type='SNV',
    #         column_specs=[
    #             # Link parent depths
    #             {'column': 'parent_DP', 'param_type': 'min',
    #              'linked_columns': ['father_DP', 'mother_DP']},
    #             # Child depth separate
    #             {'column': 'child_DP', 'param_type': 'min'},
    #             # Quality score
    #             {'column': 'quality_score', 'param_type': 'min', 'bounds': (30, 100)},
    #             # Intercept weight
    #             {'column': 'intercept', 'param_type': 'weight'},
    #         ],
    #         n_trials=500
    #     )
    # }
    # optimal_params = main(
    #     data_path='data_snv_dnms_for_filtering.tsv',
    #     reference_path='reference_parages.txt',
    #     configs=custom_configs,
    #     min_dnm_count=5,
    #     max_dnm_count=300,
    #     random_seed=RANDOM_SEED
    # )
