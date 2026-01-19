"""Plotting functions for optimization results."""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from pathlib import Path
from typing import Dict, Any, Optional, List
import logging
import optuna

logger = logging.getLogger(__name__)


def plot_optimization_results(
    data_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    best_params: Dict[str, Dict[str, Any]],
    output_dir: Optional[Path] = None,
    save_filtered: bool = True, 
    config: Optional[Any] = None
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
        filtered_data = apply_filters_from_params(data_df, var_type, best_params_var, config)

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

            # Reference data (deCODE) - include ALL samples with valid ages,
            # even those with 0 variants for this type (to match target calculation)
            all_ref_ages = reference_df[['SAMPLE', age_col]].drop_duplicates(subset=['SAMPLE']).set_index('SAMPLE')
            all_ref_ages = all_ref_ages.dropna(subset=[age_col])
            
            # Count reference variants per sample (only samples with variants will appear)
            if len(ref_subset) > 0:
                ref_counts = ref_subset.groupby('SAMPLE').size().rename('dnm_count')
            else:
                ref_counts = pd.Series(dtype=float, name='dnm_count')
            
            # Join: start with ALL samples (to include those with 0 variants for this type)
            plot_data_ref = all_ref_ages.join(ref_counts, how='left')
            plot_data_ref['dnm_count'] = plot_data_ref['dnm_count'].fillna(0)

            # Filtered input data - include ALL samples with valid ages,
            # even those with 0 variants for this type after filtering (to match optimization)
            # Get all unique samples from the FULL input data (across all variant types)
            # Use subset=['SAMPLE'] to ensure one row per sample (matching optimization)
            all_sample_ages = data_df[['SAMPLE', age_col]].drop_duplicates(subset=['SAMPLE']).set_index('SAMPLE')
            all_sample_ages = all_sample_ages.dropna(subset=[age_col])
            
            # Count filtered variants per sample (only samples with variants will appear)
            filt_counts = filtered_data.groupby('SAMPLE').size().rename('dnm_count')
            
            # Join: start with ALL samples (to include those with 0 variants for this type)
            plot_data_filt = all_sample_ages.join(filt_counts, how='left')
            plot_data_filt['dnm_count'] = plot_data_filt['dnm_count'].fillna(0)
            
            # Log sample counts for debugging
            n_with_variants = (plot_data_filt['dnm_count'] > 0).sum()
            n_total = len(plot_data_filt)
            n_zero = (plot_data_filt['dnm_count'] == 0).sum()
            mean_count = plot_data_filt['dnm_count'].mean()
            logger.info(f"  {var_type} {age_label}: {n_with_variants}/{n_total} samples have >=1 variant, {n_zero} have 0, mean={mean_count:.2f}")
            
            # Log the regression results for comparison with optimization
            if len(plot_data_filt) > 0:
                try:
                    debug_model = smf.ols(f'dnm_count ~ {age_col}', data=plot_data_filt).fit()
                    logger.info(f"  {var_type} {age_label} PLOT REGRESSION: intercept={debug_model.params['Intercept']:.4f}, slope={debug_model.params[age_col]:.4f}")
                except Exception as e:
                    logger.warning(f"  Failed to fit debug model: {e}")

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
    params: Dict[str, Any],
    config: Optional[Any] = None
) -> pd.DataFrame:
    """
    Apply filtering parameters to dataframe, handling linked columns.
    """
    # Start with variant type filter
    filtered = df[df['var_type'] == var_type].copy()

    # Build lookups for linked columns and variant_types restrictions
    linked_columns = {}
    column_variant_types = {}  # NEW: track which columns apply to which variant types
    
    if config and hasattr(config, 'optimisation'):
        for col_conf in config.optimisation.columns:
            if col_conf.linked_to:
                linked_columns[col_conf.name] = col_conf.linked_to
            if col_conf.variant_types:  # NEW: store variant type restrictions
                column_variant_types[col_conf.name] = col_conf.variant_types

    # Apply each parameter
    for param, value in params.items():
        prefix = None
        col = None
        
        # Parse the parameter name
        if param.startswith('min_'):
            prefix = 'min'
            col = param[4:]
        elif param.startswith('max_'):
            prefix = 'max'
            col = param[4:]
            
        if col:
            # Skip if this column doesn't apply to this variant type
            if col in column_variant_types:
                if var_type not in column_variant_types[col]:
                    continue  # Skip this filter for this variant type
            
            # 1. Apply to the primary column
            if col in filtered.columns:
                # Ensure column is numeric before comparison
                if not pd.api.types.is_numeric_dtype(filtered[col]):
                    filtered[col] = pd.to_numeric(filtered[col], errors='coerce')
                
                if prefix == 'min':
                    filtered = filtered[filtered[col] >= value]
                elif prefix == 'max':
                    filtered = filtered[filtered[col] <= value]
            
            # 2. Apply to the linked column (if any)
            if col in linked_columns:
                linked_col = linked_columns[col]
                if linked_col in filtered.columns:
                    # Ensure linked column is numeric before comparison
                    if not pd.api.types.is_numeric_dtype(filtered[linked_col]):
                        filtered[linked_col] = pd.to_numeric(filtered[linked_col], errors='coerce')
                    
                    if prefix == 'min':
                        filtered = filtered[filtered[linked_col] >= value]
                    elif prefix == 'max':
                        filtered = filtered[filtered[linked_col] <= value]
    
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


def plot_optuna_diagnostics(
    studies: Dict[str, optuna.Study],
    output_dir: Path,
) -> None:
    """
    Generate Optuna diagnostic plots for each variant type using matplotlib.
    
    Creates:
    - Optimization history (loss over trials)
    - Parameter importance (based on correlation with objective)
    - Parameter distributions for best trials
    - Slice plots for each parameter
    
    Parameters
    ----------
    studies : dict
        Dictionary mapping variant type to Optuna study
    output_dir : Path
        Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for var_type, study in studies.items():
        if study is None or len(study.trials) == 0:
            logger.warning(f"No trials for {var_type}, skipping Optuna plots")
            continue
        
        var_dir = output_dir / f"optuna_{var_type}"
        var_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Generating Optuna diagnostics for {var_type}...")
        
        # Get trials dataframe
        try:
            trials_df = study.trials_dataframe()
            trials_df = trials_df[trials_df['state'] == 'COMPLETE'].copy()
            if len(trials_df) == 0:
                logger.warning(f"No completed trials for {var_type}")
                continue
        except Exception as e:
            logger.warning(f"Failed to get trials dataframe for {var_type}: {e}")
            continue
        
        # Get parameter columns
        param_cols = [c for c in trials_df.columns if c.startswith('params_')]
        
        # 1. Optimization history
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            ax.scatter(trials_df['number'], trials_df['value'], 
                      alpha=0.5, s=20, c='steelblue', label='Trials')
            
            # Running minimum
            running_min = trials_df['value'].expanding().min()
            ax.plot(trials_df['number'], running_min, 
                   color='red', linewidth=2, label='Best so far')
            
            ax.set_xlabel('Trial Number', fontsize=12)
            ax.set_ylabel('Objective Value (Loss)', fontsize=12)
            ax.set_title(f'{var_type}: Optimization History', fontsize=14, fontweight='bold')
            ax.set_yscale('log')
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(var_dir / "optimization_history.png", dpi=150, bbox_inches='tight')
            plt.close()
            logger.debug(f"  Saved optimization_history.png")
        except Exception as e:
            logger.warning(f"  Failed to plot optimization history: {e}")
        
        # 2. Parameter importance (correlation-based)
        if len(param_cols) > 0:
            try:
                fig, ax = plt.subplots(figsize=(10, max(6, len(param_cols) * 0.5)))
                
                # Calculate correlation of each parameter with objective
                importances = {}
                for col in param_cols:
                    param_name = col.replace('params_', '')
                    # Use absolute correlation as importance
                    corr = trials_df[col].corr(trials_df['value'])
                    if not np.isnan(corr):
                        importances[param_name] = abs(corr)
                
                if importances:
                    # Sort by importance
                    sorted_params = sorted(importances.items(), key=lambda x: x[1], reverse=True)
                    names, values = zip(*sorted_params)
                    
                    y_pos = np.arange(len(names))
                    ax.barh(y_pos, values, color='steelblue', alpha=0.8)
                    ax.set_yticks(y_pos)
                    ax.set_yticklabels(names)
                    ax.set_xlabel('Importance (|correlation| with objective)', fontsize=12)
                    ax.set_title(f'{var_type}: Parameter Importance', fontsize=14, fontweight='bold')
                    ax.grid(True, alpha=0.3, axis='x')
                    
                    plt.tight_layout()
                    plt.savefig(var_dir / "param_importances.png", dpi=150, bbox_inches='tight')
                    plt.close()
                    logger.debug(f"  Saved param_importances.png")
            except Exception as e:
                logger.warning(f"  Failed to plot param importances: {e}")
        
        # 3. Slice plots (parameter value vs objective)
        if len(param_cols) > 0:
            try:
                n_params = len(param_cols)
                n_cols = min(3, n_params)
                n_rows = (n_params + n_cols - 1) // n_cols
                
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
                if n_params == 1:
                    axes = np.array([axes])
                axes = axes.flatten()
                
                for idx, col in enumerate(param_cols):
                    ax = axes[idx]
                    param_name = col.replace('params_', '')
                    
                    # Color by trial number
                    scatter = ax.scatter(trials_df[col], trials_df['value'],
                                        c=trials_df['number'], cmap='viridis',
                                        alpha=0.6, s=20)
                    
                    # Mark best trial
                    best_idx = trials_df['value'].idxmin()
                    ax.scatter(trials_df.loc[best_idx, col], trials_df.loc[best_idx, 'value'],
                              color='red', s=100, marker='*', zorder=5, label='Best')
                    
                    ax.set_xlabel(param_name, fontsize=10)
                    ax.set_ylabel('Objective', fontsize=10)
                    ax.set_title(param_name, fontsize=11, fontweight='bold')
                    ax.set_yscale('log')
                    ax.grid(True, alpha=0.3)
                
                # Hide unused axes
                for idx in range(n_params, len(axes)):
                    axes[idx].set_visible(False)
                
                fig.suptitle(f'{var_type}: Parameter Slice Plots', fontsize=14, fontweight='bold')
                plt.tight_layout()
                plt.savefig(var_dir / "slice_plots.png", dpi=150, bbox_inches='tight')
                plt.close()
                logger.debug(f"  Saved slice_plots.png")
            except Exception as e:
                logger.warning(f"  Failed to plot slice: {e}")
        
        # 4. Parallel coordinates plot
        if len(param_cols) >= 2:
            try:
                fig, ax = plt.subplots(figsize=(12, 6))
                
                # Normalize parameters to [0, 1] for visualization
                norm_df = trials_df[param_cols].copy()
                for col in param_cols:
                    col_min, col_max = norm_df[col].min(), norm_df[col].max()
                    if col_max > col_min:
                        norm_df[col] = (norm_df[col] - col_min) / (col_max - col_min)
                    else:
                        norm_df[col] = 0.5
                
                # Color by objective value (log scale)
                log_values = np.log10(trials_df['value'].clip(lower=1e-10))
                norm_values = (log_values - log_values.min()) / (log_values.max() - log_values.min() + 1e-10)
                
                # Plot lines
                cmap = plt.cm.viridis_r
                for idx in range(len(norm_df)):
                    color = cmap(norm_values.iloc[idx])
                    ax.plot(range(len(param_cols)), norm_df.iloc[idx].values,
                           color=color, alpha=0.3, linewidth=1)
                
                # Highlight best trial
                best_idx = trials_df['value'].idxmin()
                best_row = norm_df.loc[best_idx]
                ax.plot(range(len(param_cols)), best_row.values,
                       color='red', linewidth=3, label='Best trial')
                
                ax.set_xticks(range(len(param_cols)))
                ax.set_xticklabels([c.replace('params_', '') for c in param_cols], 
                                  rotation=45, ha='right')
                ax.set_ylabel('Normalized Parameter Value', fontsize=12)
                ax.set_title(f'{var_type}: Parallel Coordinates', fontsize=14, fontweight='bold')
                ax.legend(loc='upper right')
                ax.grid(True, alpha=0.3)
                
                # Add colorbar
                sm = plt.cm.ScalarMappable(cmap=cmap, 
                                          norm=plt.Normalize(vmin=log_values.min(), vmax=log_values.max()))
                sm.set_array([])
                cbar = plt.colorbar(sm, ax=ax, label='log10(Objective)')
                
                plt.tight_layout()
                plt.savefig(var_dir / "parallel_coordinate.png", dpi=150, bbox_inches='tight')
                plt.close()
                logger.debug(f"  Saved parallel_coordinate.png")
            except Exception as e:
                logger.warning(f"  Failed to plot parallel coordinate: {e}")
        
        # 5. Pairwise parameter scatter (contour alternative)
        if len(param_cols) >= 2:
            try:
                # Select top 4 most important parameters
                importances = {}
                for col in param_cols:
                    corr = trials_df[col].corr(trials_df['value'])
                    if not np.isnan(corr):
                        importances[col] = abs(corr)
                
                top_params = sorted(importances.items(), key=lambda x: x[1], reverse=True)[:4]
                top_param_cols = [p[0] for p in top_params]
                
                if len(top_param_cols) >= 2:
                    n_params = len(top_param_cols)
                    fig, axes = plt.subplots(n_params, n_params, figsize=(3 * n_params, 3 * n_params))
                    
                    for i, col_i in enumerate(top_param_cols):
                        for j, col_j in enumerate(top_param_cols):
                            ax = axes[i, j]
                            
                            if i == j:
                                # Diagonal: histogram
                                ax.hist(trials_df[col_i], bins=20, color='steelblue', alpha=0.7)
                                ax.set_xlabel(col_i.replace('params_', ''))
                            elif i > j:
                                # Lower triangle: scatter colored by objective
                                scatter = ax.scatter(trials_df[col_j], trials_df[col_i],
                                                    c=np.log10(trials_df['value'].clip(lower=1e-10)),
                                                    cmap='viridis_r', alpha=0.5, s=15)
                                
                                # Mark best
                                best_idx = trials_df['value'].idxmin()
                                ax.scatter(trials_df.loc[best_idx, col_j], 
                                          trials_df.loc[best_idx, col_i],
                                          color='red', s=80, marker='*', zorder=5)
                                
                                ax.set_xlabel(col_j.replace('params_', ''))
                                ax.set_ylabel(col_i.replace('params_', ''))
                            else:
                                # Upper triangle: hide
                                ax.set_visible(False)
                    
                    fig.suptitle(f'{var_type}: Parameter Pairwise Scatter', fontsize=14, fontweight='bold')
                    plt.tight_layout()
                    plt.savefig(var_dir / "pairwise_scatter.png", dpi=150, bbox_inches='tight')
                    plt.close()
                    logger.debug(f"  Saved pairwise_scatter.png")
            except Exception as e:
                logger.warning(f"  Failed to plot pairwise scatter: {e}")
    
    logger.info(f"Optuna diagnostics saved to {output_dir}")


def plot_regression_trajectory(
    trial_history: Dict[str, List[Dict[str, Any]]],
    targets: Dict[str, Dict[str, float]],
    output_dir: Path,
) -> None:
    """
    Plot slope and intercept trajectory over optimization trials.
    
    Creates a 2-row grid:
    - Row 1: Slope vs trial number, Intercept vs trial number, Loss vs trial number
    - Row 2: Slope vs Intercept (Pareto frontier), Slope error vs Intercept error
    
    Parameters
    ----------
    trial_history : dict
        Dictionary mapping variant type to list of trial results.
        Each trial result should have: 'slope', 'intercept', 'loss', 'trial_number'
    targets : dict
        Dictionary mapping variant type to target values.
        Should have: 'slope', 'intercept'
    output_dir : Path
        Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    var_types = list(trial_history.keys())
    n_var_types = len(var_types)
    
    if n_var_types == 0:
        logger.warning("No trial history to plot")
        return
    
    # Create figure with subplots for each variant type
    fig, axes = plt.subplots(n_var_types, 5, figsize=(25, 5 * n_var_types))
    if n_var_types == 1:
        axes = axes.reshape(1, -1)
    
    fig.suptitle('Optimization Trajectory: Slope & Intercept', fontsize=16, fontweight='bold')
    
    for row_idx, var_type in enumerate(var_types):
        history = trial_history[var_type]
        target = targets.get(var_type, {})
        
        if not history:
            continue
        
        # Extract data
        df = pd.DataFrame(history)
        
        target_slope = target.get('slope', np.nan)
        target_intercept = target.get('intercept', np.nan)
        
        # Calculate errors
        if not np.isnan(target_slope) and target_slope != 0:
            df['slope_error'] = ((df['slope'] - target_slope) / target_slope) ** 2
        else:
            df['slope_error'] = np.nan
            
        if not np.isnan(target_intercept) and target_intercept != 0:
            df['intercept_error'] = ((df['intercept'] - target_intercept) / target_intercept) ** 2
        else:
            df['intercept_error'] = np.nan
        
        # Plot 1: Slope over trials
        ax1 = axes[row_idx, 0]
        ax1.scatter(df['trial_number'], df['slope'], alpha=0.5, s=20, c='steelblue')
        ax1.axhline(y=target_slope, color='red', linestyle='--', linewidth=2, label=f'Target: {target_slope:.3f}')
        
        # Highlight best trial
        best_idx = df['loss'].idxmin()
        ax1.scatter(df.loc[best_idx, 'trial_number'], df.loc[best_idx, 'slope'], 
                   color='green', s=100, marker='*', zorder=5, label=f'Best: {df.loc[best_idx, "slope"]:.3f}')
        
        ax1.set_xlabel('Trial Number')
        ax1.set_ylabel('Slope')
        ax1.set_title(f'{var_type}: Slope Trajectory')
        ax1.legend(loc='best', fontsize=8)
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Intercept over trials
        ax2 = axes[row_idx, 1]
        ax2.scatter(df['trial_number'], df['intercept'], alpha=0.5, s=20, c='darkorange')
        ax2.axhline(y=target_intercept, color='red', linestyle='--', linewidth=2, label=f'Target: {target_intercept:.3f}')
        ax2.scatter(df.loc[best_idx, 'trial_number'], df.loc[best_idx, 'intercept'], 
                   color='green', s=100, marker='*', zorder=5, label=f'Best: {df.loc[best_idx, "intercept"]:.3f}')
        
        ax2.set_xlabel('Trial Number')
        ax2.set_ylabel('Intercept')
        ax2.set_title(f'{var_type}: Intercept Trajectory')
        ax2.legend(loc='best', fontsize=8)
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Loss over trials (with running minimum)
        ax3 = axes[row_idx, 2]
        ax3.scatter(df['trial_number'], df['loss'], alpha=0.5, s=20, c='gray', label='Trial loss')
        
        # Running minimum
        running_min = df['loss'].expanding().min()
        ax3.plot(df['trial_number'], running_min, color='green', linewidth=2, label='Best so far')
        
        ax3.set_xlabel('Trial Number')
        ax3.set_ylabel('Loss')
        ax3.set_title(f'{var_type}: Loss Trajectory')
        ax3.set_yscale('log')
        ax3.legend(loc='best', fontsize=8)
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Slope vs Intercept (parameter space)
        ax4 = axes[row_idx, 3]
        scatter = ax4.scatter(df['slope'], df['intercept'], c=df['loss'], 
                             cmap='viridis_r', alpha=0.6, s=20)
        ax4.scatter(target_slope, target_intercept, color='red', s=150, marker='X', 
                   zorder=5, label='Target', edgecolors='white', linewidths=2)
        ax4.scatter(df.loc[best_idx, 'slope'], df.loc[best_idx, 'intercept'], 
                   color='lime', s=150, marker='*', zorder=5, label='Best', edgecolors='black', linewidths=1)
        
        plt.colorbar(scatter, ax=ax4, label='Loss')
        ax4.set_xlabel('Slope')
        ax4.set_ylabel('Intercept')
        ax4.set_title(f'{var_type}: Parameter Space')
        ax4.legend(loc='best', fontsize=8)
        ax4.grid(True, alpha=0.3)
        
        # Plot 5: Slope error vs Intercept error
        ax5 = axes[row_idx, 4]
        valid_mask = ~(df['slope_error'].isna() | df['intercept_error'].isna())
        if valid_mask.sum() > 0:
            scatter = ax5.scatter(df.loc[valid_mask, 'slope_error'], 
                                 df.loc[valid_mask, 'intercept_error'],
                                 c=df.loc[valid_mask, 'trial_number'], 
                                 cmap='plasma', alpha=0.6, s=20)
            plt.colorbar(scatter, ax=ax5, label='Trial #')
            
            # Mark best point
            ax5.scatter(df.loc[best_idx, 'slope_error'], df.loc[best_idx, 'intercept_error'],
                       color='lime', s=150, marker='*', zorder=5, label='Best', 
                       edgecolors='black', linewidths=1)
            
            # Draw iso-loss contours (slope_error + intercept_error = constant)
            max_err = max(df.loc[valid_mask, 'slope_error'].max(), 
                         df.loc[valid_mask, 'intercept_error'].max())
            for iso_loss in [0.01, 0.1, 0.5, 1.0]:
                if iso_loss < max_err:
                    x_line = np.linspace(0, iso_loss, 100)
                    y_line = iso_loss - x_line
                    ax5.plot(x_line, y_line, 'k--', alpha=0.3, linewidth=1)
        
        ax5.set_xlabel('Slope Error (normalized)')
        ax5.set_ylabel('Intercept Error (normalized)')
        ax5.set_title(f'{var_type}: Error Decomposition')
        ax5.legend(loc='best', fontsize=8)
        ax5.grid(True, alpha=0.3)
        ax5.set_xlim(left=0)
        ax5.set_ylim(bottom=0)
    
    plt.tight_layout()
    
    plot_path = output_dir / 'regression_trajectory.png'
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved regression trajectory plot to {plot_path}")


def plot_optimization_summary(
    studies: Dict[str, optuna.Study],
    trial_history: Dict[str, List[Dict[str, Any]]],
    targets: Dict[str, Dict[str, float]],
    output_dir: Path,
) -> None:
    """
    Create a comprehensive summary plot combining key metrics.
    
    Parameters
    ----------
    studies : dict
        Dictionary mapping variant type to Optuna study
    trial_history : dict
        Dictionary mapping variant type to list of trial results
    targets : dict
        Dictionary mapping variant type to target values
    output_dir : Path
        Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    var_types = list(studies.keys())
    n_var_types = len(var_types)
    
    if n_var_types == 0:
        return
    
    fig, axes = plt.subplots(2, n_var_types, figsize=(7 * n_var_types, 12))
    if n_var_types == 1:
        axes = axes.reshape(-1, 1)
    
    fig.suptitle('Optimization Summary', fontsize=18, fontweight='bold')
    
    colors = {'SNV': 'steelblue', 'Insertion': 'darkorange', 'Deletion': 'forestgreen'}
    
    for col_idx, var_type in enumerate(var_types):
        study = studies.get(var_type)
        history = trial_history.get(var_type, [])
        target = targets.get(var_type, {})
        color = colors.get(var_type, 'gray')
        
        # Top row: Convergence plot
        ax_top = axes[0, col_idx]
        
        if study and len(study.trials) > 0:
            trials_df = study.trials_dataframe()
            trials_df = trials_df.sort_values('number')
            
            # Plot all trial values
            ax_top.scatter(trials_df['number'], trials_df['value'], 
                          alpha=0.4, s=15, c=color, label='Trials')
            
            # Running minimum
            running_min = trials_df['value'].expanding().min()
            ax_top.plot(trials_df['number'], running_min, 
                       color='red', linewidth=2, label='Best so far')
            
            ax_top.set_xlabel('Trial Number', fontsize=12)
            ax_top.set_ylabel('Loss', fontsize=12)
            ax_top.set_title(f'{var_type}: Convergence', fontsize=14, fontweight='bold')
            ax_top.set_yscale('log')
            ax_top.legend(loc='upper right', fontsize=10)
            ax_top.grid(True, alpha=0.3)
            
            # Add text with final values
            best_loss = study.best_value
            ax_top.text(0.02, 0.02, f'Best Loss: {best_loss:.2e}',
                       transform=ax_top.transAxes, fontsize=10,
                       verticalalignment='bottom', fontweight='bold',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Bottom row: Achieved vs Target
        ax_bot = axes[1, col_idx]
        
        if history:
            df = pd.DataFrame(history)
            best_idx = df['loss'].idxmin()
            
            target_slope = target.get('slope', np.nan)
            target_intercept = target.get('intercept', np.nan)
            achieved_slope = df.loc[best_idx, 'slope']
            achieved_intercept = df.loc[best_idx, 'intercept']
            
            # Bar chart comparing target vs achieved
            x = np.arange(2)
            width = 0.35
            
            target_vals = [target_slope, target_intercept]
            achieved_vals = [achieved_slope, achieved_intercept]
            
            bars1 = ax_bot.bar(x - width/2, target_vals, width, label='Target', color='royalblue', alpha=0.8)
            bars2 = ax_bot.bar(x + width/2, achieved_vals, width, label='Achieved', color=color, alpha=0.8)
            
            ax_bot.set_xticks(x)
            ax_bot.set_xticklabels(['Slope', 'Intercept'], fontsize=12)
            ax_bot.set_ylabel('Value', fontsize=12)
            ax_bot.set_title(f'{var_type}: Target vs Achieved', fontsize=14, fontweight='bold')
            ax_bot.legend(loc='best', fontsize=10)
            ax_bot.grid(True, alpha=0.3, axis='y')
            
            # Add value labels on bars
            for bar, val in zip(bars1, target_vals):
                ax_bot.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                           f'{val:.3f}', ha='center', va='bottom', fontsize=9)
            for bar, val in zip(bars2, achieved_vals):
                ax_bot.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                           f'{val:.3f}', ha='center', va='bottom', fontsize=9)
            
            # Calculate and display errors
            if target_slope != 0:
                slope_err = abs(achieved_slope - target_slope) / abs(target_slope) * 100
            else:
                slope_err = np.nan
            if target_intercept != 0:
                int_err = abs(achieved_intercept - target_intercept) / abs(target_intercept) * 100
            else:
                int_err = np.nan
            
            error_text = f'Slope Error: {slope_err:.1f}%\nIntercept Error: {int_err:.1f}%'
            ax_bot.text(0.98, 0.98, error_text,
                       transform=ax_bot.transAxes, fontsize=10,
                       verticalalignment='top', horizontalalignment='right',
                       bbox=dict(boxstyle='round', facecolor='lightcoral' if max(slope_err, int_err) > 20 else 'lightgreen', alpha=0.5))
    
    plt.tight_layout()
    
    plot_path = output_dir / 'optimization_summary.png'
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved optimization summary to {plot_path}")
