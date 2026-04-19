import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import math

# ----------- CONFIGURATION -----------
MINUTES_PER_MONTH = 43200  # 30 days * 24h * 60min

# Update metric labels to match your analysis_over_time.csv
METRIC_LABELS = {
    'Total_Cells': 'Total Cells',
    'Alive_PTEN_Deleted': 'PTEN Null Cells',
    'Alive_PTEN_Normal': 'PTEN Normal Cells',
    'PTEN_Ratio': 'PTEN Null/PTEN Normal Ratio',
    'Testosterone_Depletion_Radius': 'Testosterone Depletion Radius',
    'C_avg_all': 'Clustering Index (All Cells)',
    'C_avg_PTEN_normal': 'Clustering Index (PTEN Normal)',
    'C_avg_PTEN_deleted': 'Clustering Index (PTEN Deleted)',
    'Front_Radius_All': 'Front Radius (All Cells)',
    'Front_Radius_PTEN_Normal': 'Front Radius (PTEN Normal)',
    'Front_Radius_PTEN_Deleted': 'Front Radius (PTEN Deleted)',
    'Front_Speed_All': 'Front Speed (All Cells)',
    'Front_Speed_PTEN_Normal': 'Front Speed (PTEN Normal)',
    'Front_Speed_PTEN_Deleted': 'Front Speed (PTEN Deleted)'
}

def load_masterlist(masterlist_path):
    """Load masterlist CSV file"""
    return pd.read_csv(masterlist_path)

def filter_masterlist(masterlist_df, constant_conditions):
    """Filter masterlist by constant conditions"""
    df = masterlist_df.copy()
    for key, value in constant_conditions.items():
        if key in df.columns:
            # Convert both to string for comparison to handle numeric vs string mismatches
            df[key] = df[key].astype(str)
            df = df[df[key] == str(value)]
        else:
            print(f"Warning: Column '{key}' not found in masterlist. Available columns: {df.columns.tolist()}")
    return df

def get_scanning_groups(masterlist_df, scanning_variable):
    """Group runs by scanning variable(s) and return dict: {composition_label: [run_ids]}"""
    # Ensure scanning_variable is a list
    if isinstance(scanning_variable, str):
        scanning_variable = [scanning_variable]
    
    # Check if scanning variables exist in the dataframe
    for var in scanning_variable:
        if var not in masterlist_df.columns:
            print(f"Warning: Scanning variable '{var}' not found in filtered masterlist.")
            print(f"Available columns: {masterlist_df.columns.tolist()}")
            return {}
    
    # Group by scanning variable(s)
    grouped = masterlist_df.groupby(scanning_variable)
    groups = {}
    for name, group in grouped:
        # Convert name to a tuple if it's not already
        if not isinstance(name, tuple):
            name = (name,)
        
        # Create a label for this group
        label = ','.join([f"{var}={val}" for var, val in zip(scanning_variable, name)])
        groups[label] = group['Run_ID'].astype(str).tolist()
    
    return groups

def load_run_data(run_base_dir, run_id, analysis_rel_path, max_time_min=None):
    """Load analysis data for a single run"""
    path = os.path.join(run_base_dir, f'PCa_ABM_{run_id}', analysis_rel_path)
    if not os.path.exists(path):
        print(f"Missing: {path}")
        return None
    df = pd.read_csv(path)
    if max_time_min is not None:
        df = df[df['time_min'] <= max_time_min].reset_index(drop=True)
    return df

def aggregate_runs(run_base_dir, run_ids, analysis_rel_path, metric, max_time_min):
    """Aggregate data across multiple runs"""
    dfs = []
    times_list = []
    for run_id in run_ids:
        df = load_run_data(run_base_dir, run_id, analysis_rel_path, max_time_min)
        if df is not None:
            # Compute derived metrics if needed
            if metric == 'PTEN_Ratio':
                df['PTEN_Ratio'] = df['Alive_PTEN_Deleted'] / (df['Alive_PTEN_Normal'] + 1e-8)
            dfs.append(df)
            times_list.append(df['time_min'].values)
    
    if not dfs:
        return None
    
    # Just use time points from first run (assuming all identical)
    times = dfs[0]['time_min'].values
    # Stack data directly without interpolation
    data = np.array([df[metric].values for df in dfs])
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)
    return times, mean, std, data

def plot_temporal(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp):
    """Plot temporal profile (line plot with error bands)"""
    plt.figure(figsize=(10, 6))
    
    if not groups:
        print(f"No groups to plot for {metric}")
        return
        
    for label, run_ids in groups.items():
        result = aggregate_runs(run_base_dir, run_ids, analysis_rel_path, metric, final_time_months * MINUTES_PER_MONTH)
        if result is None:
            continue
        times, mean, std, data = result
        plt.plot(times / MINUTES_PER_MONTH, mean, label=label)
        plt.fill_between(times / MINUTES_PER_MONTH, mean - std, mean + std, alpha=0.2)
    
    plt.xlabel('Time (months)')
    plt.ylabel(METRIC_LABELS.get(metric, metric))
    plt.title(f'{METRIC_LABELS.get(metric, metric)} vs Time')
    
    # Add legend only if we plotted something
    if len(plt.gca().get_lines()) > 0:
        plt.legend()
    
    plt.tight_layout()
    fname = f"{metric}__{constant_conditions}_temporal_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname))
    plt.show()

def plot_violin(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp):
    """Plot distribution at final time point (violin plot)"""
    plt.figure(figsize=(10, 6))
    
    if not groups:
        print(f"No groups to plot for {metric}")
        return
    
    violin_data = []
    violin_labels = []
    
    for label, run_ids in groups.items():
        result = aggregate_runs(run_base_dir, run_ids, analysis_rel_path, metric, final_time_months * MINUTES_PER_MONTH)
        if result is None:
            continue
        times, mean, std, data = result
        idx = np.argmin(np.abs(times / MINUTES_PER_MONTH - final_time_months))
        values = data[:, idx]
        if len(values) > 0:
            violin_data.append(values)
            violin_labels.append(label)
    
    if not violin_data:
        print(f"No data to plot violin for {metric}")
        return
        
    sns.violinplot(data=violin_data)
    plt.xticks(ticks=range(len(violin_labels)), labels=violin_labels, rotation=45, ha='right')
    plt.xlabel('Tissue Composition')
    plt.ylabel(METRIC_LABELS.get(metric, metric))
    plt.title(f'{METRIC_LABELS.get(metric, metric)} at {final_time_months} months')
    plt.tight_layout()
    fname = f"{metric}_{constant_conditions}_violin_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname))
    plt.show()

def plot_metrics(
    masterlist_path,
    run_base_dir,
    analysis_rel_path,
    constant_conditions,
    scanning_variable,
    metrics,
    plot_type='temporal',
    final_time_months=15,
    output_dir='datanalaysis_output'
):
    """Main function to generate plots based on specified conditions"""
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    print(f"Loading masterlist from: {masterlist_path}")
    masterlist_df = load_masterlist(masterlist_path)
    
    print(f"Filtering by constant conditions: {constant_conditions}")
    filtered_df = filter_masterlist(masterlist_df, constant_conditions)
    
    if filtered_df.empty:
        print("No runs match the specified conditions!")
        return
    
    print(f"Found {len(filtered_df)} matching runs")
    print(f"Grouping by scanning variables: {scanning_variable}")
    groups = get_scanning_groups(filtered_df, scanning_variable)
    
    if not groups:
        print("No valid groups found after scanning!")
        return
    
    print(f"Found {len(groups)} groups: {list(groups.keys())}")
    
    for metric in metrics:
        print(f"Plotting {metric} ({plot_type})")
        if plot_type == 'temporal':
            plot_temporal(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp)
        elif plot_type == 'violin':
            plot_violin(groups, run_base_dir, analysis_rel_path, metric, constant_conditions, final_time_months, output_dir, timestamp)
        else:
            raise ValueError("plot_type must be 'temporal' or 'violin'")
        # --- Statistical analysis ---
        print(f"Statistical analysis for {metric} at {final_time_months} months:")
        statistical_analysis_groups(groups, run_base_dir, analysis_rel_path, metric, final_time_months)


def statistical_analysis_groups(groups, run_base_dir, analysis_rel_path, metric, final_time_months):
    """
    For each group, collect metric values at the final timepoint and perform statistical tests.
    """
    from scipy.stats import ttest_ind, f_oneway, levene, bartlett
    data_dict = {}
    for label, run_ids in groups.items():
        result = aggregate_runs(run_base_dir, run_ids, analysis_rel_path, metric, final_time_months * MINUTES_PER_MONTH)
        if result is None:
            continue
        times, mean, std, data = result
        idx = np.argmin(np.abs(times / MINUTES_PER_MONTH - final_time_months))
        values = data[:, idx]
        data_dict[label] = values

    if len(data_dict) < 2:
        print(f"Not enough groups for statistical comparison for {metric}.")
        return

    # Compare means
    values = list(data_dict.values())
    labels = list(data_dict.keys())
    if len(values) == 2:
        stat, p = ttest_ind(values[0], values[1])
        print(f"T-test for {metric}: stat={stat:.3f}, p={p:.3e} ({labels[0]} vs {labels[1]})")
    else:
        stat, p = f_oneway(*values)
        print(f"ANOVA for {metric}: stat={stat:.3f}, p={p:.3e} (groups: {labels})")

    # Compare variances
    stat_lev, p_lev = levene(*values)
    stat_bart, p_bart = bartlett(*values)
    print(f"Levene's test for {metric}: stat={stat_lev:.3f}, p={p_lev:.3e}")
    print(f"Bartlett's test for {metric}: stat={stat_bart:.3f}, p={p_bart:.3e}")

def calculate_segregation_index(df):
    """
    Calculate Cohen's κ (segregation index) from clustering data.
    
    Uses the correct formula:
    κ = (p_o - p_e) / (1 - p_e)
    
    Where:
    - p_o = observed same-type edge fraction
    - p_e = expected same-type edge fraction (random mixing)
    
    p_o = f_S * C_avg_S + f_R * C_avg_R
    p_e = f_S^2 + f_R^2
    
    Where:
    - f_S = fraction of S cells (PTEN_normal)
    - f_R = fraction of R cells (PTEN_deleted)
    - C_avg_S = average clustering for S cells
    - C_avg_R = average clustering for R cells
    
    κ ranges from:
    +1 → perfect clustering (all same-type neighbors)
     0 → random mixing
    -1 → perfect alternation / checkerboard
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns: time_min, C_avg_PTEN_normal, C_avg_PTEN_deleted,
                                Alive_PTEN_Normal, Alive_PTEN_Deleted
    
    Returns:
    --------
    df : pd.DataFrame
        Input df with added columns: freq_S, freq_R, p_observed, p_expected, kappa
    """
    df = df.copy()
    
    # Calculate population frequencies
    # S = PTEN_normal (Sensitive)
    # R = PTEN_deleted (Resistant)
    total_alive = df['Alive_PTEN_Normal'] + df['Alive_PTEN_Deleted']
    df['freq_S'] = df['Alive_PTEN_Normal'] / (total_alive + 1e-10)
    df['freq_R'] = df['Alive_PTEN_Deleted'] / (total_alive + 1e-10)
    
    # Rename for clarity
    df['f_S'] = df['freq_S']
    df['f_R'] = df['freq_R']
    
    # Observed same-type edge fraction
    # p_o = f_S * C_avg_S + f_R * C_avg_R
    df['p_observed'] = (df['f_S'] * df['C_avg_PTEN_normal'] + 
                        df['f_R'] * df['C_avg_PTEN_deleted'])
    
    # Expected same-type edge fraction (random mixing)
    # p_e = f_S^2 + f_R^2
    df['p_expected'] = df['f_S']**2 + df['f_R']**2
    
    # Cohen's κ (segregation index)
    # κ = (p_o - p_e) / (1 - p_e)
    # Avoid division by zero when p_e = 1 (only one type present)
    df['kappa'] = np.where(
        df['p_expected'] < 1.0,
        (df['p_observed'] - df['p_expected']) / (1.0 - df['p_expected']),
        0.0
    )
    
    return df

def postprocess_segregation_indices(groups, run_base_dir, analysis_rel_path, max_time_min, output_dir, timestamp):
    """
    Post-process all runs to calculate segregation indices and save to new CSV files.
    
    Parameters:
    -----------
    groups : dict
        Dictionary mapping group labels to list of run_ids
    run_base_dir : str
        Base directory containing PCa_ABM_* directories
    analysis_rel_path : str
        Relative path to analysis_over_time.csv within each run directory
    max_time_min : float
        Maximum simulation time in minutes to include
    output_dir : str
        Output directory for processed files
    timestamp : str
        Timestamp string for naming output files
    """
    os.makedirs(output_dir, exist_ok=True)
    
    for label, run_ids in groups.items():
        print(f"\nProcessing segregation indices for group: {label}")
        
        group_dfs = []
        for run_id in run_ids:
            df = load_run_data(run_base_dir, run_id, analysis_rel_path, max_time_min)
            if df is not None:
                # Calculate segregation indices
                df = calculate_segregation_index(df)
                df['Run_ID'] = run_id
                group_dfs.append(df)
                print(f"  Run {run_id}: processed successfully")
            else:
                print(f"  Run {run_id}: failed to load")
        
        if group_dfs:
            # Concatenate all runs in this group
            combined_df = pd.concat(group_dfs, ignore_index=True)
            
            # Save to CSV
            safe_label = label.replace('=', '_').replace(',', '_')
            output_file = os.path.join(output_dir, f"segregation_indices_{safe_label}_{timestamp}.csv")
            combined_df.to_csv(output_file, index=False)
            print(f"  Saved to: {output_file}")
            
            # Print summary statistics
            print(f"\n  Summary statistics for {label}:")
            print(f"    κ (mean): {combined_df['kappa'].mean():.3f} ± {combined_df['kappa'].std():.3f}")
            print(f"    κ (min): {combined_df['kappa'].min():.3f}, κ (max): {combined_df['kappa'].max():.3f}")
            print(f"    p_observed (mean): {combined_df['p_observed'].mean():.3f}")
            print(f"    p_expected (mean): {combined_df['p_expected'].mean():.3f}")

def plot_segregation_index_temporal(groups, run_base_dir, analysis_rel_path, max_time_min, output_dir, timestamp):
    """
    Plot temporal evolution of segregation indices using Cohen's κ formula.
    """
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    for label, run_ids in groups.items():
        dfs = []
        for run_id in run_ids:
            df = load_run_data(run_base_dir, run_id, analysis_rel_path, max_time_min)
            if df is not None:
                df = calculate_segregation_index(df)
                dfs.append(df)
        
        if dfs:
            times = dfs[0]['time_min'].values
            
            # Cohen's κ
            kappa_data = np.array([df['kappa'].values for df in dfs])
            mean_kappa = np.mean(kappa_data, axis=0)
            std_kappa = np.std(kappa_data, axis=0)
            
            axes[0].plot(times / MINUTES_PER_MONTH, mean_kappa, label=label, linewidth=2)
            axes[0].fill_between(times / MINUTES_PER_MONTH, 
                           mean_kappa - std_kappa, 
                           mean_kappa + std_kappa, 
                           alpha=0.2)
            
            # p_observed vs p_expected
            p_obs_data = np.array([df['p_observed'].values for df in dfs])
            p_exp_data = np.array([df['p_expected'].values for df in dfs])
            mean_p_obs = np.mean(p_obs_data, axis=0)
            mean_p_exp = np.mean(p_exp_data, axis=0)
            
            axes[1].plot(times / MINUTES_PER_MONTH, mean_p_obs, label=f'{label} (observed)', linewidth=2)
            axes[1].plot(times / MINUTES_PER_MONTH, mean_p_exp, label=f'{label} (expected)', 
                        linewidth=2, linestyle='--')
    
    # κ plot
    axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axes[0].set_xlabel('Time (months)', fontsize=12)
    axes[0].set_ylabel("Cohen's κ (Segregation Index)", fontsize=12)
    axes[0].set_title("Cohen's κ: Temporal Evolution", fontsize=14)
    axes[0].legend(fontsize=10)
    axes[0].grid(alpha=0.3)
    
    # p_observed vs p_expected plot
    axes[1].set_xlabel('Time (months)', fontsize=12)
    axes[1].set_ylabel("Same-type edge fraction", fontsize=12)
    axes[1].set_title("Observed vs Expected Same-type Edge Fraction", fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    
    fname = f"segregation_index_temporal_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname), dpi=300)
    print(f"Saved figure: {fname}")
    plt.show()

def plot_segregation_violin(groups, run_base_dir, analysis_rel_path, final_time_months, output_dir, timestamp):
    """
    Create violin plots for segregation metrics at final timepoint.
    
    Plots:
    1. Cohen's κ (segregation index)
    2. p_observed (observed same-type edge fraction)
    3. p_expected (expected same-type edge fraction by random mixing)
    """
    metrics_to_plot = ['kappa', 'p_observed', 'p_expected']
    metric_labels = {
        'kappa': "Cohen's κ (Segregation Index)",
        'p_observed': 'Observed Same-type Edge Fraction',
        'p_expected': 'Expected Same-type Edge Fraction'
    }
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    for idx, metric in enumerate(metrics_to_plot):
        violin_data = []
        violin_labels = []
        
        for label, run_ids in groups.items():
            group_values = []
            
            for run_id in run_ids:
                df = load_run_data(run_base_dir, run_id, analysis_rel_path, final_time_months * MINUTES_PER_MONTH)
                if df is not None:
                    # Calculate segregation indices
                    df = calculate_segregation_index(df)
                    
                    # Get value at final timepoint
                    final_value = df[metric].iloc[-1]
                    group_values.append(final_value)
            
            if len(group_values) > 0:
                violin_data.append(group_values)
                violin_labels.append(label)
        
        if violin_data:
            # Create violin plot
            parts = axes[idx].violinplot(violin_data, positions=range(len(violin_data)), 
                                         showmeans=True, showmedians=True)
            
            axes[idx].set_xticks(range(len(violin_labels)))
            axes[idx].set_xticklabels(violin_labels, rotation=45, ha='right')
            axes[idx].set_ylabel(metric_labels[metric], fontsize=11)
            axes[idx].set_title(metric_labels[metric], fontsize=12, fontweight='bold')
            axes[idx].grid(axis='y', alpha=0.3)
            
            # Add reference lines
            if metric == 'kappa':
                axes[idx].axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=2, label='Random mixing')
                axes[idx].legend()
            elif metric in ['p_observed', 'p_expected']:
                axes[idx].set_ylim([0, 1])
    
    plt.tight_layout()
    
    fname = f"segregation_violin_final_timepoint_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches='tight')
    print(f"Saved figure: {fname}")
    plt.show()


# def plot_segregation_violin_seaborn(groups, run_base_dir, analysis_rel_path, final_time_months, output_dir, timestamp):
#     """
#     Create violin plots using seaborn for better aesthetics.
    
#     Creates a combined dataframe and plots all three metrics.
#     Groups are ordered by cell adhesion multiplier strength.
#     """
#     # Collect all data
#     all_data = []
    
#     for label, run_ids in groups.items():
#         for run_id in run_ids:
#             df = load_run_data(run_base_dir, run_id, analysis_rel_path, final_time_months * MINUTES_PER_MONTH)
#             if df is not None:
#                 # Calculate segregation indices
#                 df = calculate_segregation_index(df)
                
#                 # Get final timepoint
#                 final_row = df.iloc[-1]
                
#                 all_data.append({
#                     'Group': label,
#                     'Run_ID': run_id,
#                     "Cohen's κ": final_row['kappa'],
#                     'p_observed': final_row['p_observed'],
#                     'p_expected': final_row['p_expected'],
#                     'f_S': final_row['f_S'],
#                     'f_R': final_row['f_R']
#                 })
    
#     if not all_data:
#         print("No data to plot!")
#         return
    
#     data_df = pd.DataFrame(all_data)
    
#     # Extract adhesion multiplier values and sort
#     # Expected format: "Cell_cell_adhesion_multiplier=X"
#     adhesion_order = [0.1, 0.5, 1, 3, 7, 10]
    
#     # Create mapping of group labels to adhesion values
#     group_adhesion_map = {}
#     for group in data_df['Group'].unique():
#         try:
#             adhesion_val = float(group.split('=')[-1])
#             group_adhesion_map[group] = adhesion_val
#         except:
#             group_adhesion_map[group] = float('inf')  # Put unparseable groups at end
    
#     # Sort data_df by adhesion strength
#     data_df['adhesion_value'] = data_df['Group'].map(group_adhesion_map)
#     data_df = data_df.sort_values('adhesion_value')
    
#     # Create ordered group list for plotting
#     ordered_groups = data_df['Group'].unique()
    
#     # Create subplots
#     fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
#     # Plot 1: Cohen's κ
#     sns.violinplot(data=data_df, x='Group', y="Cohen's κ", ax=axes[0, 0], palette='Set2', order=ordered_groups)
#     axes[0, 0].axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Random mixing')
#     axes[0, 0].set_ylabel("Cohen's κ", fontsize=11)
#     axes[0, 0].set_title("Cohen's κ (Segregation Index)", fontsize=12, fontweight='bold')
#     axes[0, 0].tick_params(axis='x', rotation=45)
#     axes[0, 0].legend()
#     axes[0, 0].grid(axis='y', alpha=0.3)
    
#     # Plot 2: p_observed
#     sns.violinplot(data=data_df, x='Group', y='p_observed', ax=axes[0, 1], palette='Set2', order=ordered_groups)
#     axes[0, 1].set_ylabel('Same-type Edge Fraction', fontsize=11)
#     axes[0, 1].set_title('Observed Same-type Edge Fraction (p_o)', fontsize=12, fontweight='bold')
#     axes[0, 1].set_ylim([0, 1])
#     axes[0, 1].tick_params(axis='x', rotation=45)
#     axes[0, 1].grid(axis='y', alpha=0.3)
    
#     # Plot 3: p_expected
#     sns.violinplot(data=data_df, x='Group', y='p_expected', ax=axes[1, 0], palette='Set2', order=ordered_groups)
#     axes[1, 0].set_ylabel('Same-type Edge Fraction', fontsize=11)
#     axes[1, 0].set_title('Expected Same-type Edge Fraction by Random Mixing (p_e)', fontsize=12, fontweight='bold')
#     axes[1, 0].set_ylim([0, 1])
#     axes[1, 0].tick_params(axis='x', rotation=45)
#     axes[1, 0].grid(axis='y', alpha=0.3)
    
#     # Plot 4: p_observed vs p_expected comparison
#     comparison_data = []
#     for _, row in data_df.iterrows():
#         comparison_data.append({
#             'Group': row['Group'], 
#             'p': row['p_observed'], 
#             'Type': 'Observed',
#             'adhesion_value': row['adhesion_value']
#         })
#         comparison_data.append({
#             'Group': row['Group'], 
#             'p': row['p_expected'], 
#             'Type': 'Expected',
#             'adhesion_value': row['adhesion_value']
#         })
    
#     comparison_df = pd.DataFrame(comparison_data)
#     sns.violinplot(data=comparison_df, x='Group', y='p', hue='Type', ax=axes[1, 1], 
#                    split=True, palette=['skyblue', 'lightcoral'], order=ordered_groups)
#     axes[1, 1].set_ylabel('Same-type Edge Fraction', fontsize=11)
#     axes[1, 1].set_title('Observed vs Expected: Same-type Edge Fraction', fontsize=12, fontweight='bold')
#     axes[1, 1].set_ylim([0, 1])
#     axes[1, 1].tick_params(axis='x', rotation=45)
#     axes[1, 1].grid(axis='y', alpha=0.3)
#     axes[1, 1].legend(title='Type', fontsize=10)
    
#     plt.tight_layout()
    
#     fname = f"segregation_violin_seaborn_{timestamp}.png"
#     plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches='tight')
#     print(f"Saved figure: {fname}")
#     plt.show()
    
#     # Print summary statistics in order of adhesion strength
#     print("\n" + "="*60)
#     print("SEGREGATION METRICS AT FINAL TIMEPOINT")
#     print("="*60)
#     for group in ordered_groups:
#         group_data = data_df[data_df['Group'] == group]
#         kappa_col = "Cohen's κ"
#         print(f"\n{group}:")
#         print(f"  Cohen's κ: {group_data[kappa_col].mean():.3f} ± {group_data[kappa_col].std():.3f}")
#         print(f"  p_observed: {group_data['p_observed'].mean():.3f} ± {group_data['p_observed'].std():.3f}")
#         print(f"  p_expected: {group_data['p_expected'].mean():.3f} ± {group_data['p_expected'].std():.3f}")
    
#     return data_df
def plot_segregation_violin_seaborn(groups, run_base_dir, analysis_rel_path, final_time_months, output_dir, timestamp):
    """
    Create violin plots using seaborn for better aesthetics.
    
    Creates a combined dataframe and plots all three metrics.
    Groups are ordered by cell adhesion multiplier strength.
    """
    # Collect all data
    all_data = []
    
    for label, run_ids in groups.items():
        for run_id in run_ids:
            df = load_run_data(run_base_dir, run_id, analysis_rel_path, final_time_months * MINUTES_PER_MONTH)
            if df is not None:
                # Calculate segregation indices
                df = calculate_segregation_index(df)
                
                # Get final timepoint
                final_row = df.iloc[-1]
                
                all_data.append({
                    'Group': label,
                    'Run_ID': run_id,
                    "Cohen's κ": final_row['kappa'],
                    'p_observed': final_row['p_observed'],
                    'p_expected': final_row['p_expected'],
                    'f_S': final_row['f_S'],
                    'f_R': final_row['f_R']
                })
    
    if not all_data:
        print("No data to plot!")
        return
    
    data_df = pd.DataFrame(all_data)
    
    # Extract adhesion multiplier values and sort
    # Expected format: "Cell_cell_adhesion_multiplier=X"
    adhesion_order = [0.1, 0.5, 1, 3, 7, 10]
    
    # Create mapping of group labels to adhesion values
    group_adhesion_map = {}
    for group in data_df['Group'].unique():
        try:
            adhesion_val = float(group.split('=')[-1])
            group_adhesion_map[group] = adhesion_val
        except:
            group_adhesion_map[group] = float('inf')  # Put unparseable groups at end
    
    # Sort data_df by adhesion strength
    data_df['adhesion_value'] = data_df['Group'].map(group_adhesion_map)
    data_df = data_df.sort_values('adhesion_value')
    
    # Create ordered group list for plotting
    ordered_groups = data_df['Group'].unique()
    
    # Create subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Cohen's κ
    sns.violinplot(data=data_df, x='Group', y="Cohen's κ", ax=axes[0, 0], palette='Set2', order=ordered_groups)
    axes[0, 0].axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Random mixing')
    axes[0, 0].set_ylabel("Cohen's κ", fontsize=11)
    axes[0, 0].set_title("Cohen's κ (Segregation Index)", fontsize=12, fontweight='bold')
    axes[0, 0].set_xticklabels(ordered_groups, rotation=45, ha='right')
    axes[0, 0].legend()
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    # Plot 2: p_observed
    sns.violinplot(data=data_df, x='Group', y='p_observed', ax=axes[0, 1], palette='Set2', order=ordered_groups)
    axes[0, 1].set_ylabel('Same-type Edge Fraction', fontsize=11)
    axes[0, 1].set_title('Observed Same-type Edge Fraction (p_o)', fontsize=12, fontweight='bold')
    axes[0, 1].set_xticklabels(ordered_groups, rotation=45, ha='right')
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    # Plot 3: p_expected
    sns.violinplot(data=data_df, x='Group', y='p_expected', ax=axes[1, 0], palette='Set2', order=ordered_groups)
    axes[1, 0].set_ylabel('Same-type Edge Fraction', fontsize=11)
    axes[1, 0].set_title('Expected Same-type Edge Fraction by Random Mixing (p_e)', fontsize=12, fontweight='bold')
    axes[1, 0].set_xticklabels(ordered_groups, rotation=45, ha='right')
    axes[1, 0].grid(axis='y', alpha=0.3)
    
    # Plot 4: p_observed vs p_expected comparison
    comparison_data = []
    for _, row in data_df.iterrows():
        comparison_data.append({
            'Group': row['Group'], 
            'p': row['p_observed'], 
            'Type': 'Observed',
            'adhesion_value': row['adhesion_value']
        })
        comparison_data.append({
            'Group': row['Group'], 
            'p': row['p_expected'], 
            'Type': 'Expected',
            'adhesion_value': row['adhesion_value']
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    sns.violinplot(data=comparison_df, x='Group', y='p', hue='Type', ax=axes[1, 1], 
                   split=True, palette=['skyblue', 'lightcoral'], order=ordered_groups)
    axes[1, 1].set_ylabel('Same-type Edge Fraction', fontsize=11)
    axes[1, 1].set_title('Observed vs Expected: Same-type Edge Fraction', fontsize=12, fontweight='bold')
    axes[1, 1].set_xticklabels(ordered_groups, rotation=45, ha='right')
    axes[1, 1].grid(axis='y', alpha=0.3)
    axes[1, 1].legend(title='Type', fontsize=10)
    
    plt.tight_layout()
    
    fname = f"segregation_violin_seaborn_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches='tight')
    print(f"Saved figure: {fname}")
    plt.show()
    
    # Print summary statistics in order of adhesion strength
    print("\n" + "="*60)
    print("SEGREGATION METRICS AT FINAL TIMEPOINT")
    print("="*60)
    for group in ordered_groups:
        group_data = data_df[data_df['Group'] == group]
        kappa_col = "Cohen's κ"
        print(f"\n{group}:")
        print(f"  Cohen's κ: {group_data[kappa_col].mean():.3f} ± {group_data[kappa_col].std():.3f}")
        print(f"  p_observed: {group_data['p_observed'].mean():.3f} ± {group_data['p_observed'].std():.3f}")
        print(f"  p_expected: {group_data['p_expected'].mean():.3f} ± {group_data['p_expected'].std():.3f}")
    
    return data_df

def plot_clustering_vs_adhesion(
    masterlist_path,
    run_base_dir,
    analysis_rel_path,
    cohorts,
    androgen_condition='High',
    scenario='Scenario0',
    uptake_rate_multiplier=1,
    final_time_months=15,
    output_dir='datanalaysis_output'
):
    """
    Plot clustering index (C_avg_PTEN_deleted) vs log10 adhesion strength
    for multiple cohorts.
    
    Parameters:
    -----------
    masterlist_path : str
        Path to masterlist CSV
    run_base_dir : str
        Base directory containing PCa_ABM_* directories
    analysis_rel_path : str
        Relative path to analysis_over_time.csv
    cohorts : list
        List of cohorts to plot (e.g., ['BR', 'CTRL', 'TR'])
    androgen_condition : str
        Androgen condition to filter (default: 'High')
    scenario : str
        Scenario to filter (default: 'Scenario0')
    uptake_rate_multiplier : float
        Uptake rate multiplier to filter (default: 1)
    final_time_months : float
        Final timepoint to extract (default: 15)
    output_dir : str
        Output directory for plots
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    print(f"Loading masterlist from: {masterlist_path}")
    masterlist_df = load_masterlist(masterlist_path)
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Color palette for cohorts
    colors = {'BR': '#1f77b4', 'CTRL': '#ff7f0e', 'TR': '#2ca02c'}
    markers = {'BR': 'o', 'CTRL': 's', 'TR': '^'}
    
    for cohort in cohorts:
        print(f"\nProcessing cohort: {cohort}")
        
        # Filter by cohort, androgen condition, scenario, and uptake rate multiplier
        filtered_df = filter_masterlist(
            masterlist_df,
            {
                'Cohort': cohort,
                'Androgen_condition': androgen_condition,
                'Scenario': scenario,
                'Uptake_rate_multiplier': uptake_rate_multiplier
            }
        )
        
        if filtered_df.empty:
            print(f"  No runs found for {cohort}")
            continue
        
        print(f"  Found {len(filtered_df)} runs for {cohort}")
        print(f"  Columns: {filtered_df.columns.tolist()}")
        
        # Get unique adhesion multipliers
        if 'Cell_cell_adhesion_multiplier' not in filtered_df.columns:
            print(f"  'Cell_cell_adhesion_multiplier' column not found")
            print(f"  Available columns: {filtered_df.columns.tolist()}")
            continue
        
        # Convert to numeric and filter out non-numeric values
        adhesion_values_raw = filtered_df['Cell_cell_adhesion_multiplier'].unique()
        print(f"  Raw adhesion values: {adhesion_values_raw}")
        
        adhesion_values = []
        for val in adhesion_values_raw:
            try:
                num_val = float(val)
                adhesion_values.append(num_val)
            except (ValueError, TypeError):
                print(f"    Skipping non-numeric value: {val}")
                pass
        
        if not adhesion_values:
            print(f"  No valid numeric adhesion values found. Available values were: {adhesion_values_raw}")
            continue
        
        adhesion_values = sorted(adhesion_values)
        print(f"  Found {len(adhesion_values)} numeric adhesion values: {adhesion_values}")
        
        clustering_means = []
        clustering_stds = []
        log_adhesion_values = []
        
        for adhesion in adhesion_values:
            # Filter runs for this adhesion value - need to convert column to float for comparison
            try:
                filtered_df_numeric = filtered_df.copy()
                filtered_df_numeric['Cell_cell_adhesion_multiplier'] = pd.to_numeric(filtered_df_numeric['Cell_cell_adhesion_multiplier'], errors='coerce')
                subset_df = filtered_df_numeric[filtered_df_numeric['Cell_cell_adhesion_multiplier'] == adhesion]
            except Exception as e:
                print(f"    Adhesion {adhesion}: Error converting to numeric - {e}")
                continue
            
            if subset_df.empty:
                print(f"    Adhesion {adhesion}: No matching runs")
                continue
            
            run_ids = subset_df['Run_ID'].astype(str).tolist()
            print(f"    Adhesion {adhesion}: Found {len(run_ids)} runs")
            
            # Aggregate clustering data
            result = aggregate_runs(
                run_base_dir,
                run_ids,
                analysis_rel_path,
                'C_avg_PTEN_deleted',
                final_time_months * MINUTES_PER_MONTH
            )
            
            if result is None:
                print(f"    Adhesion {adhesion}: No data could be aggregated")
                continue
            
            times, mean, std, data = result
            
            # Get value at final timepoint
            final_idx = len(mean) - 1
            final_value = mean[final_idx]
            final_std = std[final_idx]
            
            clustering_means.append(final_value)
            clustering_stds.append(final_std)
            log_adhesion_values.append(math.log10(adhesion))
            
            print(f"    Adhesion {adhesion} (log10={math.log10(adhesion):.2f}): C_avg_PTEN_deleted = {final_value:.3f} ± {final_std:.3f}")
        
        if clustering_means:
            # Plot with error bars
            ax.errorbar(
                log_adhesion_values,
                clustering_means,
                yerr=clustering_stds,
                marker=markers[cohort],
                linestyle='-',
                linewidth=2,
                markersize=8,
                label=cohort,
                color=colors[cohort],
                capsize=5,
                alpha=0.8
            )
    
    # Formatting
    ax.set_xlabel('log10(Cell-cell Adhesion Multiplier)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Clustering Index (C_avg_PTEN_deleted)', fontsize=13, fontweight='bold')
    ax.set_title(f'Clustering Index vs Adhesion Strength\n({androgen_condition} AR, {scenario})',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    fname = f"clustering_vs_adhesion_{timestamp}.png"
    plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches='tight')
    print(f"\nSaved figure: {os.path.join(output_dir, fname)}")
    plt.show()

# ----------- EXAMPLE USAGE -----------
if __name__ == '__main__':
    # Plot clustering index vs adhesion strength for all three cohorts
    # Filtered for Uptake_rate_multiplier = 1
    plot_clustering_vs_adhesion(
        masterlist_path='/jet/home/sskemkar/PCa_code/ABMruns_masterlist_prostatecancer.csv',
        run_base_dir='/ocean/projects/mcb200052p/sskemkar',
        analysis_rel_path='output/analysis_over_time.csv',
        cohorts=['BR', 'CTRL', 'TR'],
        androgen_condition='High',
        scenario='Scenario2',
        uptake_rate_multiplier=1,
        final_time_months=15,
        output_dir='mar2026'
    )

    # NEW: Post-process segregation indices
    # print("\n" + "="*60)
    # print("POST-PROCESSING: CALCULATING SEGREGATION INDICES")
    # print("="*60)
    
    # masterlist_df = load_masterlist('/jet/home/sskemkar/PCa_code/ABMruns_masterlist_prostatecancer.csv')
    # filtered_df = filter_masterlist(masterlist_df, {'Cohort': 'CTRL', 'Scenario': 'Scenario2', 'Androgen_condition': 'High', 'Uptake_rate_multiplier': '1', 'Uptake_celltype': 'both'})
    # groups = get_scanning_groups(filtered_df, ['Cell_cell_adhesion_multiplier'])
    
    # timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # # Calculate and save segregation indices
    # postprocess_segregation_indices(
    #     groups=groups,
    #     run_base_dir='/ocean/projects/mcb200052p/sskemkar',
    #     analysis_rel_path='output/analysis_over_time.csv',
    #     max_time_min=15 * MINUTES_PER_MONTH,
    #     output_dir='datanalaysis_output_oct212025',
    #     timestamp=timestamp
    # )
    
    # # Plot segregation index temporal evolution
    # plot_segregation_index_temporal(
    #     groups=groups,
    #     run_base_dir='/ocean/projects/mcb200052p/sskemkar',
    #     analysis_rel_path='output/analysis_over_time.csv',
    #     max_time_min=15 * MINUTES_PER_MONTH,
    #     output_dir='datanalaysis_output_oct212025',
    #     timestamp=timestamp
    # )

    # # Violin plots at final timepoint (RECOMMENDED - better aesthetics)
    # plot_segregation_violin_seaborn(
    #     groups=groups,
    #     run_base_dir='/ocean/projects/mcb200052p/sskemkar',
    #     analysis_rel_path='output/analysis_over_time.csv',
    #     final_time_months=15,
    #     output_dir='datanalaysis_output_oct212025',
    #     timestamp=timestamp
    # )
    
    # # Alternative: simpler violin plots
    # plot_segregation_violin(
    #     groups=groups,
    #     run_base_dir='/ocean/projects/mcb200052p/sskemkar',
    #     analysis_rel_path='output/analysis_over_time.csv',
    #     final_time_months=15,
    #     output_dir='datanalaysis_output_oct212025',
    #     timestamp=timestamp
    # )
