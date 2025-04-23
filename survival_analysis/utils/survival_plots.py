import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import logging
import sys


def plot_forest(df, show_proteins="all", sort_by="hazard_ratio", plot_title="Cox Proportional Hazard Analysis", significance_column="FDR_adjusted_pval",save_path=None):
    """
    Create a Cox proportional hazard forest plot.

    Parameters:
    - df: DataFrame containing columns ['protein', 'hazard_ratio', 'CI_lower', 'CI_upper', 'FDR_adjusted_pval', 'p_value']
    - show_proteins: "all" to display all proteins, or a list of specific proteins to filter before plotting.
    - sort_by: Column name to sort proteins by ("hazard_ratio", "p_value", or "FDR_adjusted_pval").
    - plot_title: Custom title for the plot.
    - significance_column: Column to use for significance annotation ("FDR_adjusted_pval" or "p_value").

    Returns:
    - A forest plot displaying hazard ratios and confidence intervals.
    """

    # Validate sorting input
    valid_sorting_options = ["hazard_ratio", "p_value", "FDR_adjusted_pval"]
    if sort_by not in valid_sorting_options:
        raise ValueError(f"Invalid sort_by option. Choose from {valid_sorting_options}")

    # Validate significance column input
    if significance_column not in ["FDR_adjusted_pval", "p_value"]:
        raise ValueError(f"Invalid significance_column option. Choose between 'FDR_adjusted_pval' or 'p_value'")

    # Filter dataset based on selected proteins
    if show_proteins != "all":
        df = df[df["protein"].isin(show_proteins)]

    # Sorting logic
    if sort_by == "hazard_ratio":
        df_sorted = pd.concat([
            df[(df["hazard_ratio"] > 1) & (df["protein"] != "LRG1")].sort_values(by="hazard_ratio", ascending=False),
            df[df["protein"] == "LRG1"],
            df[df["hazard_ratio"] < 1].sort_values(by="hazard_ratio", ascending=True)
        ])
    else:
        df_sorted = df.sort_values(by=sort_by, ascending=True)  # Sort by chosen metric

    # Create the forest plot
    plt.figure(figsize=(20, len(df_sorted) * 0.8))

    # Plot error bars (confidence intervals)
    plt.hlines(y=np.arange(len(df_sorted)), xmin=df_sorted["CI_lower"], xmax=df_sorted["CI_upper"], color="black")

    # Scatter plot for hazard ratios with colors based on conditions
    for i, (prot, hr) in enumerate(zip(df_sorted["protein"], df_sorted["hazard_ratio"])):
        
        plt.scatter(hr, i, color='black', edgecolor='black', s=20, zorder=4, marker='o')

    # Add reference line at HR = 1
    plt.axvline(x=1, color='steelblue', linestyle='--', linewidth=1, alpha=0.7)  

    # Convert adjusted p-values to asterisks for significance levels
    def get_significance(pval):
        if pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return ''

    # Annotate significance levels above the dot, using the selected column
    for i, (pval, hr) in enumerate(zip(df_sorted[significance_column], df_sorted["hazard_ratio"])):
        significance = get_significance(pval)
        plt.text(hr, i + 0.1, significance, horizontalalignment='center', fontsize=10, color='black', fontweight='bold')

   
    # Remove spines for a cleaner look
    sns.despine(top=True, right=True, left=False, bottom=False)

    # Formatting
    plt.yticks(np.arange(len(df_sorted)), df_sorted["protein"], fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlabel('Hazard Ratio (95% CI)\n Survival from Enrollment', fontsize=12)
    plt.title(f'{plot_title}', fontsize=14)
    plt.grid(False)

    # Add a legend for significance levels
    legend_text = "* p < 0.05  |  ** p < 0.01  |  *** p < 0.001"
    plt.text(2/9, 0.4, legend_text, fontsize=12, color='black', verticalalignment='top')

    # Set x-axis range
    plt.xlim(0.2, 2)

    # Save
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)  # Create directory if it doesn't exist
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    # ✅ Show plot after saving (prevents blank saved images)
    plt.show()
    plt.close() 


def plot_forest_two_cohorts_with_statistics(discovery_df, replication_df, proteins_to_plot, save_path=None):
    """
    Create a Cox Proportional Hazard forest plot for selected proteins from discovery & replication cohorts.

    Parameters:
    - discovery_df: DataFrame containing CoxPH results for the discovery cohort.
    - replication_df: DataFrame containing CoxPH results for the replication cohort.
    - proteins_to_plot: List of protein names to visualize.
    - save_path: Optional, file path to save the figure. If None, the plot is displayed without saving.

    Returns:
    - A professional forest plot comparing hazard ratios across discovery and replication cohorts.
    """

    # Filter for selected proteins
    df_dis = discovery_df[discovery_df["protein"].isin(proteins_to_plot)].copy()
    df_rep = replication_df[replication_df["protein"].isin(proteins_to_plot)].copy()

    # Assign cohort labels
    df_dis["cohort"] = df_dis["protein"] + " (Discovery)"
    df_rep["cohort"] = df_rep["protein"] + " (Replication)"

    # Merge datasets for visualization
    df_combined = pd.concat([df_dis, df_rep])

    # Extract values for plotting
    hazard_ratios = df_combined["hazard_ratio"]
    proteins = df_combined["cohort"]
    ci_lower = df_combined["CI_lower"]
    ci_upper = df_combined["CI_upper"]
    concordance_index = df_combined["concordance_index"]

    # **Determine the correct significance column**
    significance_values = []
    for index, row in df_combined.iterrows():
        if "Discovery" in row["cohort"]:
            significance_values.append(row["FDR_adjusted_pval"])  # Use FDR for Discovery
        else:
            significance_values.append(row["p_value"])  # Use raw p-value for Replication

    # **Plot Setup**
    plt.figure(figsize=(16, len(proteins_to_plot) * 3))

    # Plot error bars (confidence intervals)
    plt.hlines(y=np.arange(len(proteins), 0, -1), xmin=ci_lower, xmax=ci_upper, color="black", linewidth=2)

    # Scatter plot for hazard ratios
    plt.scatter(hazard_ratios, np.arange(len(proteins), 0, -1), color="black", edgecolor="black", s=200, zorder=3, marker="o")

    # Add reference line at HR = 1
    plt.axvline(x=1, color='steelblue', linestyle='--', linewidth=2, alpha=0.7)

    # Function to convert p-values to significance stars
    def get_significance(pval):
        if pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return ''

    # Annotate significance levels above the dot using **correct p-value source**
    for i, (pval, hr) in enumerate(zip(significance_values, hazard_ratios)):
        significance = get_significance(pval)
        plt.text(hr, len(proteins) - i, significance, horizontalalignment='center', fontsize=22, color='black', fontweight='bold')

    # Annotate CI, adjusted p-values, and concordance index
    for i, (ci_l, ci_u, pval, cindex, hr) in enumerate(zip(ci_lower, ci_upper, significance_values, concordance_index, hazard_ratios)):
        plt.text(ci_u + 0.05, len(proteins) - i, 
                 f"HR={hr:.2f} (CI 95%: [{ci_l:.2f}, {ci_u:.2f}])\nCI-index={cindex:.3f}",
                 verticalalignment='center', fontsize=16, color="black", fontweight="bold")

    # Formatting
    sns.despine(top=True, right=True, left=False, bottom=False)
    plt.yticks(np.arange(len(proteins), 0, -1), proteins, fontsize=18)
    plt.xticks(fontsize=18)
    plt.xlabel('Hazard Ratio', fontsize=20)
    plt.title('CoxPH univariate', fontsize=30, fontweight="bold")

    # Set x-axis range
    plt.xlim(0.5, 3)

    # Improve layout and visibility
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

    # **Save or Show the Plot**
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)  # Ensure directory exists
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()
    plt.close()



#---------------------------------------------------------------------------

logging.basicConfig(level=logging.INFO)

def find_best_shared_split_and_plot_km(discovery_df, replication_df, proteins, time_col, event_col, min_samples=50):
    """
    Finds the best shared split point that maximizes log-rank test significance in both discovery & replication cohorts.
    Ensures at least `min_samples` in each group.
    Performs Kaplan-Meier analysis for Discovery, Replication, and Combined cohorts and plots KM curves.

    Parameters:
    - discovery_df: DataFrame for the discovery cohort.
    - replication_df: DataFrame for the replication cohort.
    - proteins: List of proteins to analyze.
    - time_col: Column name representing survival time.
    - event_col: Column name for event occurrence (1=event, 0=censored).
    - min_samples: Minimum number of samples required in each group (default=50).

    Returns:
    - DataFrame with split points and log-rank statistics.
    """

    results_list = []

    for protein in proteins:
        logging.info(f"Processing protein: {protein}")

        # Remove missing values
        valid_discovery = discovery_df.dropna(subset=[protein, time_col, event_col])
        valid_replication = replication_df.dropna(subset=[protein, time_col, event_col])

        # Unique values sorted
        unique_values = np.sort(valid_discovery[protein].unique())

        best_split = None
        best_statistic_combined = -np.inf

        # Iterate through potential split points
        for split_point in unique_values:
            group1_dis = valid_discovery[valid_discovery[protein] <= split_point]
            group2_dis = valid_discovery[valid_discovery[protein] > split_point]
            group1_rep = valid_replication[valid_replication[protein] <= split_point]
            group2_rep = valid_replication[valid_replication[protein] > split_point]

            # Ensure each group has at least `min_samples` in both cohorts
            if len(group1_dis) >= min_samples and len(group2_dis) >= min_samples and \
               len(group1_rep) >= min_samples and len(group2_rep) >= min_samples:

                # Log-rank tests
                logrank_dis = logrank_test(group1_dis[time_col], group2_dis[time_col],
                                           event_observed_A=group1_dis[event_col], event_observed_B=group2_dis[event_col])
                logrank_rep = logrank_test(group1_rep[time_col], group2_rep[time_col],
                                           event_observed_A=group1_rep[event_col], event_observed_B=group2_rep[event_col])

                # Combined score
                combined_statistic = logrank_dis.test_statistic + logrank_rep.test_statistic

                # Store the best split point
                if combined_statistic > best_statistic_combined and logrank_dis.p_value < 0.05 and logrank_rep.p_value < 0.05:
                    best_split = split_point
                    best_statistic_combined = combined_statistic

        if best_split is not None:
            # Apply the shared best split
            df_dis_low = valid_discovery[valid_discovery[protein] <= best_split]
            df_dis_high = valid_discovery[valid_discovery[protein] > best_split]
            df_rep_low = valid_replication[valid_replication[protein] <= best_split]
            df_rep_high = valid_replication[valid_replication[protein] > best_split]

            # Combine Discovery + Replication
            df_combined = pd.concat([valid_discovery, valid_replication])
            df_comb_low = df_combined[df_combined[protein] <= best_split]
            df_comb_high = df_combined[df_combined[protein] > best_split]

            # Fit Kaplan-Meier models
            kmf_dis_low, kmf_dis_high = KaplanMeierFitter(), KaplanMeierFitter()
            kmf_rep_low, kmf_rep_high = KaplanMeierFitter(), KaplanMeierFitter()
            kmf_comb_low, kmf_comb_high = KaplanMeierFitter(), KaplanMeierFitter()

            kmf_dis_low.fit(df_dis_low[time_col], event_observed=df_dis_low[event_col])
            kmf_dis_high.fit(df_dis_high[time_col], event_observed=df_dis_high[event_col])
            kmf_rep_low.fit(df_rep_low[time_col], event_observed=df_rep_low[event_col])
            kmf_rep_high.fit(df_rep_high[time_col], event_observed=df_rep_high[event_col])
            kmf_comb_low.fit(df_comb_low[time_col], event_observed=df_comb_low[event_col])
            kmf_comb_high.fit(df_comb_high[time_col], event_observed=df_comb_high[event_col])

            # Log-rank tests
            logrank_dis = logrank_test(df_dis_low[time_col], df_dis_high[time_col],
                                       event_observed_A=df_dis_low[event_col], event_observed_B=df_dis_high[event_col])
            logrank_rep = logrank_test(df_rep_low[time_col], df_rep_high[time_col],
                                       event_observed_A=df_rep_low[event_col], event_observed_B=df_rep_high[event_col])
            logrank_comb = logrank_test(df_comb_low[time_col], df_comb_high[time_col],
                                        event_observed_A=df_comb_low[event_col], event_observed_B=df_comb_high[event_col])

            # Store results
            results_list.append({
                'Protein': protein,
                'Best Shared Split Point': best_split,
                'Discovery Log-Rank χ²': logrank_dis.test_statistic,
                'Discovery p-value': logrank_dis.p_value,
                'Replication Log-Rank χ²': logrank_rep.test_statistic,
                'Replication p-value': logrank_rep.p_value,
                'Combined Log-Rank χ²': logrank_comb.test_statistic,
                'Combined p-value': logrank_comb.p_value,
                'Samples Discovery (Low)': len(df_dis_low),
                'Samples Discovery (High)': len(df_dis_high),
                'Samples Replication (Low)': len(df_rep_low),
                'Samples Replication (High)': len(df_rep_high),
                'Samples Combined (Low)': len(df_comb_low),
                'Samples Combined (High)': len(df_comb_high)
            })

            # Plot Kaplan-Meier curves
            plot_km_curve(kmf_dis_low, kmf_dis_high, df_dis_low, df_dis_high, protein, "Discovery Cohort", logrank_dis)
            plot_km_curve(kmf_rep_low, kmf_rep_high, df_rep_low, df_rep_high, protein, "Replication Cohort", logrank_rep)
            plot_km_curve(kmf_comb_low, kmf_comb_high, df_comb_low, df_comb_high, protein, "Combined Cohort", logrank_comb)

    return pd.DataFrame(results_list)

# Function to plot Kaplan-Meier curve
def plot_km_curve(kmf_low, kmf_high, df_low, df_high, protein_name, title, logrank):
    plt.figure(figsize=(8, 4))
    
    # Plot survival curves
    kmf_low.plot(ci_show=False, linewidth=2, color='steelblue')
    kmf_high.plot(ci_show=False, linewidth=2, color='crimson')

    plt.gca().get_legend().remove()

    # Log-rank test annotation
    plt.text(60, 0.7, f'log-rank χ²={logrank.test_statistic:.1f}, P={logrank.p_value:.4f}', fontsize=10)

    # Median survival annotation
    plt.text(68, 0.6, f'Low {protein_name}, median: {df_low["Survival_from_onset (months)"].median():.1f} m (n={len(df_low)})', 
             fontsize=10, color='steelblue')
    plt.text(68, 0.5, f'High {protein_name}, median: {df_high["Survival_from_onset (months)"].median():.1f} m (n={len(df_high)})', 
             fontsize=10, color='crimson')

    plt.xlabel("Survival from onset (months)", size=14)
    plt.ylabel("Percent survival(%)", size=14)
    plt.title(f"{title}", fontsize=14, fontweight="bold")
    sns.despine(top=True, right=True, left=False, bottom=False)

    plt.show()
