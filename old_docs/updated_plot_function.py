"""Updated plotting function using new parameter names from dnm-harmoniser."""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_results(ukbb_df, decode_df, all_best_params):
    """
    Plot DNMs vs mid-parental age using optimized parameters.

    Parameters updated to match new dnm-harmoniser column naming:
    - Uses exact column names from configuration
    - Uses max_* for minimum optimization, min_* for maximum optimization
    - Range constraints have both min_* and max_*
    """
    if ukbb_df is None or decode_df is None or not all_best_params:
        print("\nSkipping plot due to missing data or optimization results.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(24, 7), sharey=False)
    fig.suptitle('De Novo Mutations vs. Mid-Parental Age by Variant Type',
                 fontsize=20, fontweight='bold')

    var_types_to_plot = ['SNV', 'Insertion', 'Deletion']

    for i, var_type in enumerate(var_types_to_plot):
        ax = axes[i]
        best_params = all_best_params.get(var_type)

        if not best_params:
            print(f"No optimal parameters found for {var_type}, skipping plot.")
            continue

        # Print parameters for debugging
        print(f"\n{var_type} parameters:")
        for key, value in sorted(best_params.items()):
            print(f"  {key}: {value}")

        # Apply filters using NEW parameter names
        filters = [ukbb_df['var_type'] == var_type]

        # Minimum optimization columns (use max_* thresholds - filter OUT high values)
        if 'max_DeNovoCNN_prob' in best_params:
            filters.append(ukbb_df['DeNovoCNN_prob'] <= best_params['max_DeNovoCNN_prob'])

        if 'max_DNM' in best_params:
            filters.append(ukbb_df['DNM'] <= best_params['max_DNM'])

        if 'max_MQ' in best_params:
            filters.append(ukbb_df['MQ'] <= best_params['max_MQ'])

        if 'max_nparAADn0' in best_params:
            filters.append(ukbb_df['nparAADn0'] <= best_params['max_nparAADn0'])

        # Maximum optimization columns (use min_* thresholds - filter OUT low values)
        if 'min_child_coverage' in best_params:
            filters.append(ukbb_df['child_coverage'] >= best_params['min_child_coverage'])

        if 'min_father_coverage' in best_params:
            filters.append(ukbb_df['father_coverage'] >= best_params['min_father_coverage'])

        if 'min_mother_coverage' in best_params:
            filters.append(ukbb_df['mother_coverage'] >= best_params['min_mother_coverage'])

        # Range optimization columns (use both min_* and max_*)
        # VAF - symmetric range
        if 'min_VAF' in best_params and 'max_VAF' in best_params:
            filters.append(ukbb_df['VAF'] >= best_params['min_VAF'])
            filters.append(ukbb_df['VAF'] <= best_params['max_VAF'])

        # IMF - regular range (only for indels)
        if var_type in ['Insertion', 'Deletion']:
            if 'min_IMF' in best_params and 'max_IMF' in best_params:
                imf_float = ukbb_df['IMF'].astype('float')
                filters.append(imf_float >= best_params['min_IMF'])
                filters.append(imf_float <= best_params['max_IMF'])

        # Handle FS if it exists in parameters (you may need to add to config)
        if 'min_FS' in best_params:
            filters.append(ukbb_df['FS'] >= best_params['min_FS'])

        # Apply all filters
        filter_mask = pd.Series([True] * len(ukbb_df))
        for f in filters:
            filter_mask &= f

        filtered_ukbb_df = ukbb_df[filter_mask].copy()

        print(f"  Filtered {len(ukbb_df)} → {len(filtered_ukbb_df)} variants ({len(filtered_ukbb_df)/len(ukbb_df)*100:.1f}%)")

        # Prepare deCODE data for this variant type
        decode_subset = decode_df[decode_df['var_type'] == var_type]
        decode_counts = decode_subset.groupby('pid').size().rename('dnm_count')
        decode_ages = decode_subset[['pid', 'midparage']].drop_duplicates().set_index('pid')
        plot_data_decode = decode_ages.join(decode_counts, how='left').fillna(0)

        # Prepare UKBB data for this variant type
        ukbb_counts = filtered_ukbb_df.groupby('SAMPLE').size().rename('dnm_count')
        ukbb_ages = filtered_ukbb_df[['SAMPLE', 'midparage']].drop_duplicates().set_index('SAMPLE')
        plot_data_ukbb = ukbb_ages.join(ukbb_counts, how='left').fillna(0)

        # Plotting
        sns.regplot(data=plot_data_decode, x='midparage', y='dnm_count',
                    label='deCODE', color='royalblue', ax=ax,
                    scatter_kws={'alpha': 0.6, 's': 40, 'edgecolor': 'w'},
                    line_kws={'linewidth': 2})
        sns.regplot(data=plot_data_ukbb, x='midparage', y='dnm_count',
                    label='UKBB (Filtered)', color='darkorange', ax=ax,
                    scatter_kws={'alpha': 0.7, 's': 40, 'edgecolor': 'w'},
                    line_kws={'linewidth': 2})

        ax.set_title(f'{var_type}s', fontsize=16)
        ax.set_xlabel('Mid-Parental Age', fontsize=14)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend()

    axes[0].set_ylabel('Number of DNMs per Person', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


# Example usage:
if __name__ == "__main__":
    print("Example parameter dictionary structure:")
    print("""
    all_best_params = {
        'SNV': {
            'max_DeNovoCNN_prob': 0.5,
            'max_DNM': 10.0,
            'max_MQ': 40,
            'max_nparAADn0': 5,
            'min_child_coverage': 10,
            'min_father_coverage': 8,
            'min_mother_coverage': 8,
            'min_VAF': 30,
            'max_VAF': 70,
        },
        'Insertion': {
            # ... same parameters ...
            'min_IMF': 0.3,
            'max_IMF': 0.7,
        },
        'Deletion': {
            # ... same parameters ...
            'min_IMF': 0.3,
            'max_IMF': 0.7,
        }
    }
    """)
