"""Plotting functions for optimization results."""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from pathlib import Path
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


def plot_optimization_results(
    data_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    best_params: Dict[str, Dict[str, Any]],
    output_dir: Optional[Path] = None,
    save_filtered: bool = True
) -> None:
    """
    Plot DNMs vs. parental age separately for each parent and variant type.

    Creates a 2x3 grid:
    - Row 1: Paternal age vs DNMs (SNV, Insertion, Deletion)
    - Row 2: Maternal age vs DNMs (SNV, Insertion, Deletion)

    Parameters
    ----------
    data_df : pd.DataFrame
        Input data (e.g., UKBB)
    reference_df : pd.DataFrame
        Reference data (e.g., deCODE)
    best_params : dict
        Best parameters for each variant type
    output_dir : Path, optional
        Directory to save plots and filtered data
    save_filtered : bool
        Whether to save filtered variants to file
    """
    if data_df is None or reference_df is None or not best_params:
        logger.warning("Skipping plot due to missing data or optimization results.")
        return

    # Create output directory if specified
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Create 2x3 subplot grid
    fig, axes = plt.subplots(2, 3, figsize=(24, 14))
    fig.suptitle('De Novo Mutations vs. Parental Age by Variant Type',
                 fontsize=20, fontweight='bold', y=0.995)

    var_types = ['SNV', 'Insertion', 'Deletion']
    age_types = [
        ('paternal_age', 'Paternal Age', 0),
        ('maternal_age', 'Maternal Age', 1)
    ]

    # Store all filtered data
    all_filtered_variants = []

    for col_idx, var_type in enumerate(var_types):
        best_params_var = best_params.get(var_type)

        if not best_params_var:
            logger.warning(f"No optimal parameters found for {var_type}, skipping plot.")
            for row_idx in [0, 1]:
                axes[row_idx, col_idx].text(
                    0.5, 0.5, f'No parameters\nfor {var_type}',
                    ha='center', va='center', fontsize=14
                )
                axes[row_idx, col_idx].set_title(f'{var_type}s', fontsize=16)
            continue

        # Apply filters to get filtered data
        filtered_data = apply_filters_from_params(data_df, var_type, best_params_var)

        # Save filtered variants
        if save_filtered and output_dir:
            filtered_copy = filtered_data.copy()
            filtered_copy['variant_type'] = var_type
            all_filtered_variants.append(filtered_copy)

        logger.info(f"{var_type}: {len(data_df[data_df['var_type']==var_type])} → "
                   f"{len(filtered_data)} variants "
                   f"({len(filtered_data)/max(len(data_df[data_df['var_type']==var_type]),1)*100:.1f}%)")

        # Prepare reference data for this variant type
        ref_subset = reference_df[reference_df['var_type'] == var_type]

        # Plot for each age type
        for age_col, age_label, row_idx in age_types:
            ax = axes[row_idx, col_idx]

            # Reference data (deCODE) - exclude samples with NA ages
            ref_counts = ref_subset.groupby('SAMPLE').size().rename('dnm_count')
            ref_ages = ref_subset[['SAMPLE', age_col]].drop_duplicates().set_index('SAMPLE')
            plot_data_ref = ref_ages.join(ref_counts, how='left')
            plot_data_ref['dnm_count'] = plot_data_ref['dnm_count'].fillna(0)
            plot_data_ref = plot_data_ref.dropna(subset=[age_col])

            # Filtered input data - exclude samples with NA ages
            filt_counts = filtered_data.groupby('SAMPLE').size().rename('dnm_count')
            filt_ages = filtered_data[['SAMPLE', age_col]].drop_duplicates().set_index('SAMPLE')
            plot_data_filt = filt_ages.join(filt_counts, how='left')
            plot_data_filt['dnm_count'] = plot_data_filt['dnm_count'].fillna(0)
            plot_data_filt = plot_data_filt.dropna(subset=[age_col])

            # Plot scatter points only (no automatic regression lines)
            if len(plot_data_ref) > 0:
                ax.scatter(
                    plot_data_ref[age_col], plot_data_ref['dnm_count'],
                    label='Reference (deCODE)', color='royalblue',
                    alpha=0.6, s=40, edgecolor='w'
                )

            if len(plot_data_filt) > 0:
                ax.scatter(
                    plot_data_filt[age_col], plot_data_filt['dnm_count'],
                    label='Filtered Data', color='darkorange',
                    alpha=0.7, s=40, edgecolor='w'
                )

            # Fit and plot actual regression lines using the coefficients
            # Fit reference regression
            if len(plot_data_ref) > 0:
                try:
                    ref_model = smf.ols(f'dnm_count ~ {age_col}', data=plot_data_ref).fit()
                    age_range = np.linspace(plot_data_ref[age_col].min(), plot_data_ref[age_col].max(), 100)
                    ref_pred = ref_model.params['Intercept'] + ref_model.params[age_col] * age_range
                    ax.plot(age_range, ref_pred, color='royalblue', linewidth=2.5,
                           label=f'Reference fit (slope={ref_model.params[age_col]:.3f})')
                except:
                    pass

            # Fit filtered data regression
            if len(plot_data_filt) > 0:
                try:
                    filt_model = smf.ols(f'dnm_count ~ {age_col}', data=plot_data_filt).fit()
                    age_range = np.linspace(plot_data_filt[age_col].min(), plot_data_filt[age_col].max(), 100)
                    filt_pred = filt_model.params['Intercept'] + filt_model.params[age_col] * age_range
                    ax.plot(age_range, filt_pred, color='darkorange', linewidth=2.5,
                           label=f'Filtered fit (slope={filt_model.params[age_col]:.3f})')
                except:
                    pass

            # Formatting
            ax.set_title(f'{var_type}s', fontsize=16, fontweight='bold')
            ax.set_xlabel(age_label, fontsize=14)
            ax.set_ylabel('Number of DNMs per Person' if col_idx == 0 else '', fontsize=14)
            ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)
            ax.legend(loc='upper left', fontsize=11)

    plt.tight_layout()

    # Save plot if output directory specified
    if output_dir:
        plot_path = output_dir / 'optimization_results.png'
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        logger.debug(f"Saved plot to {plot_path}")

    # Close figure to free memory (don't display)
    plt.close()

    # Save filtered variants to file
    if save_filtered and output_dir and all_filtered_variants:
        all_filtered_df = pd.concat(all_filtered_variants, ignore_index=True)
        output_file = output_dir / 'filtered_variants.tsv'
        all_filtered_df.to_csv(output_file, sep='\t', index=False)
        logger.debug(f"Saved {len(all_filtered_df)} filtered variants to {output_file}")

        # Save summary statistics
        summary_file = output_dir / 'filter_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("FILTERING SUMMARY\n")
            f.write("=" * 80 + "\n\n")

            for var_type in var_types:
                if var_type in best_params:
                    original_count = len(data_df[data_df['var_type'] == var_type])
                    filtered_count = len(all_filtered_df[all_filtered_df['variant_type'] == var_type])
                    pct = filtered_count / max(original_count, 1) * 100

                    f.write(f"{var_type}:\n")
                    f.write(f"  Original: {original_count:,}\n")
                    f.write(f"  Filtered: {filtered_count:,}\n")
                    f.write(f"  Retained: {pct:.1f}%\n\n")

                    f.write(f"  Parameters:\n")
                    for param, value in sorted(best_params[var_type].items()):
                        f.write(f"    {param}: {value}\n")
                    f.write("\n")

        logger.debug(f"Saved summary to {summary_file}")


def apply_filters_from_params(
    df: pd.DataFrame,
    var_type: str,
    params: Dict[str, Any]
) -> pd.DataFrame:
    """
    Apply filtering parameters to dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    var_type : str
        Variant type to filter
    params : dict
        Filtering parameters (e.g., {'max_MQ': 40, 'min_child_coverage': 10})

    Returns
    -------
    pd.DataFrame
        Filtered dataframe
    """
    # Start with variant type filter
    filtered = df[df['var_type'] == var_type].copy()

    # Apply each parameter
    for param, value in params.items():
        if param.startswith('min_'):
            col = param[4:]  # Remove 'min_' prefix
            if col in filtered.columns:
                filtered = filtered[filtered[col] >= value]

        elif param.startswith('max_'):
            col = param[4:]  # Remove 'max_' prefix
            if col in filtered.columns:
                filtered = filtered[filtered[col] <= value]

    return filtered


def save_parameters(
    best_params: Dict[str, Dict[str, Any]],
    output_dir: Path,
    filename: str = 'best_parameters.txt'
) -> None:
    """
    Save best parameters to a human-readable file.

    Parameters
    ----------
    best_params : dict
        Best parameters for each variant type
    output_dir : Path
        Output directory
    filename : str
        Output filename
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / filename

    with open(output_file, 'w') as f:
        f.write("OPTIMAL FILTERING PARAMETERS\n")
        f.write("=" * 80 + "\n\n")

        for var_type, params in best_params.items():
            f.write(f"{var_type}:\n")
            f.write("-" * 40 + "\n")

            for param, value in sorted(params.items()):
                if isinstance(value, float):
                    f.write(f"  {param}: {value:.4f}\n")
                else:
                    f.write(f"  {param}: {value}\n")

            f.write("\n")

    logger.debug(f"Saved parameters to {output_file}")
