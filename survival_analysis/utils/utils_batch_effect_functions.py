
import os
import sys
sys.path.insert(0, '/home/labs/hornsteinlab/hadarkl/Hadar/primals/code')

import pandas as pd
import numpy as np

#from neuroCombat import neuroCombat
from pycombat import pycombat
from pycombat import Combat

import scanpy as sc
from scipy.stats import ttest_ind,shapiro, levene
from statsmodels.stats.multitest import multipletests

from scipy.stats import shapiro, levene
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler, MinMaxScaler

#import umap
from configs.constants import MAIN_DIR, seed,MODEL_DIR,DISCOVERY_INDEX_BATCH1



#plot 2 distributions dogether 
# def plot_combined_distributions(discovery_df, replication_df):
#     # Fill NaN values with a small number close to zero
#     discovery_df = discovery_df.fillna(1e-10)
#     replication_df = replication_df.fillna(1e-10)

#     # Calculate row-wise means for each sample
#     discovery_mean = discovery_df.mean(axis=1)
#     replication_mean = replication_df.mean(axis=1)

#     # Create a combined DataFrame for plotting
#     combined_data = pd.DataFrame({
#         'Value': pd.concat([discovery_mean, replication_mean], axis=0),
#         'Batch': ['Discovery'] * len(discovery_mean) + ['Replication'] * len(replication_mean)
#     })

#     # Plot overlapping distributions using seaborn
#     plt.figure(figsize=(8, 6))
#     sns.kdeplot(data=combined_data, x='Value', hue='Batch', fill=True, common_norm=False, alpha=0.5)
#     plt.title('Distribution of Discovery and Replication Samples after batch correction[combat]')
#     plt.xlabel('Mean Protein Abundance')
#     plt.ylabel('Density')
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.tight_layout()
#     plt.show()


def plot_combined_distributions(df, df_combat):
    # Calculate row-wise means for each sample
    df_means = df.mean(axis=1)  # Mean before imputation
    df_combat_means = df_combat.mean(axis=1)  # Mean after imputation

    # Create a combined DataFrame for plotting
    combined_data = pd.DataFrame({
        'Value': pd.concat([df_means, df_combat_means], axis=0),
        'Data': ['Before combat'] * len(df_means) + ['After combat'] * len(df_combat_means)
    })

    # Plot overlapping distributions using seaborn
    plt.figure(figsize=(8, 6))
    sns.kdeplot(data=combined_data, x='Value', hue='Data', fill=True, common_norm=False, alpha=0.5)
    plt.title('Distribution before and after combat')
    plt.xlabel('Mean Protein Abundance')
    plt.ylabel('Density')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

def boxplot_DE_proteins_after_batch_correction(protein_list, discovery_df, replication_df, group_size=20, colors=None, figsize=(6, 4)):
    """
    Plot the distribution of proteins in groups for discovery and replication datasets
    and calculate the medians and absolute differences.

    Parameters:
    - protein_list: List of proteins to plot.
    - discovery_df: DataFrame containing discovery data.
    - replication_df: DataFrame containing replication data.
    - group_size: Number of proteins per plot.
    - colors: Dictionary mapping 'Discovery' and 'Replication' to custom colors.
    - figsize: Tuple specifying figure size.

    Returns:
    - A tuple containing:
      - A list of median results as dictionaries.
      - A DataFrame with median values and absolute differences for each protein.
    """
    

    # Filter proteins to ensure they exist in both datasets
    valid_proteins = [protein for protein in protein_list if protein in discovery_df.columns and protein in replication_df.columns]
    if not valid_proteins:
        print("No valid proteins found in both datasets.")
        return [], pd.DataFrame()

    # Default colors if none are provided
    if colors is None:
        colors = {'Discovery': 'steelblue', 'Replication': 'lightcoral'}

    # Initialize a list to store median differences
    median_results = []

    # Split the valid proteins into groups
    for i in range(0, len(valid_proteins), group_size):
        subset_proteins = valid_proteins[i:i+group_size]

        # Prepare the data
        discovery_long = discovery_df[subset_proteins].melt(var_name='Protein', value_name='Value')
        discovery_long['Type'] = 'Discovery'

        replication_long = replication_df[subset_proteins].melt(var_name='Protein', value_name='Value')
        replication_long['Type'] = 'Replication'

        combined_df = pd.concat([discovery_long, replication_long])

        # Plot with custom colors
        plt.figure(figsize=figsize)
        sns.boxplot(data=combined_df, x='Protein', y='Value', hue='Type', palette=colors)
        plt.title(f"Protein Distribution after Combat (Proteins {i+1} to {i+len(subset_proteins)})", fontsize=16)
        plt.xlabel("Proteins", fontsize=20)
        plt.ylabel("Expression Level", fontsize=20)
        plt.legend(title="Type", fontsize=12, loc='upper right')
        plt.xticks(rotation=45, fontsize=15)
        plt.tight_layout()
        plt.show()

        # Calculate medians and absolute differences
        for protein in subset_proteins:
            discovery_median = discovery_df[protein].median()
            replication_median = replication_df[protein].median()
            median_difference = abs(discovery_median - replication_median)
            median_results.append({
                'Protein': protein,
                'Discovery Median': discovery_median,
                'Replication Median': replication_median,
                'Absolute Median Difference': median_difference
            })

    # Convert results to a DataFrame
    median_results_df = pd.DataFrame(median_results)
    return median_results, median_results_df


def plot_protein_groups(protein_list, discovery_df, replication_df, group_size=20):
    """
    Plot the distribution of proteins in groups for discovery and replication datasets.

    Parameters:
    - protein_list: List of proteins to plot.
    - discovery_df: DataFrame containing discovery data.
    - replication_df: DataFrame containing replication data.
    - group_size: Number of proteins per plot.

    Returns:
    - Generates multiple plots, each showing the distribution of a subset of proteins.
    """
    # Split the proteins into groups
    for i in range(0, len(protein_list), group_size):
        subset_proteins = protein_list[i:i+group_size]

        # Prepare the data
        discovery_long = discovery_df[subset_proteins].melt(var_name='Protein', value_name='Value')
        discovery_long['Type'] = 'Discovery'

        replication_long = replication_df[subset_proteins].melt(var_name='Protein', value_name='Value')
        replication_long['Type'] = 'Replication'

        combined_df = pd.concat([discovery_long, replication_long])

        # Plot
        plt.figure(figsize=(6, 4))
        sns.boxplot(data=combined_df, x='Protein', y='Value', hue='Type', palette="Set2")
        plt.title(f"Protein Distribution (Proteins {i+1} to {i+len(subset_proteins)})", fontsize=16)
        plt.xlabel("Proteins", fontsize=14)
        plt.ylabel("Expression Level", fontsize=14)
        plt.legend(title="Type", fontsize=12, loc='upper right')
        plt.xticks(rotation=45, fontsize=10)
        plt.tight_layout()
        plt.show()

# pca plot

def PCA_plot(df, n_components=2):
    """
    Perform PCA on the given DataFrame, print explained variance, and plot the results 
    with points colored by the 'batch' column (discovery = blue, replication = red).
    
    """
    # Separate features and batch column
    features = df.drop(columns=['Batch'])  # Exclude the 'batch' column
    batch_labels = df['Batch']
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(features)
    
    # Print explained variance
    explained_variance = pca.explained_variance_ratio_ * 100  # Convert to percentage
    for i, var in enumerate(explained_variance):
        print(f"Explained Variance by PC{i+1}: {var:.2f}%")
    
    # Create a DataFrame for the PCA results
    pca_df = pd.DataFrame(data=principal_components, columns=[f'PC{i+1}' for i in range(n_components)])
    pca_df['Batch'] = batch_labels.values  # Add the batch column
    
    # Plot the PCA results
    plt.figure(figsize=(8, 4))
    for batch, color in zip(['discovery', 'replication'], ['steelblue', 'lightcoral']):
        subset = pca_df[pca_df['Batch'] == batch]
        plt.scatter(subset['PC1'], subset['PC2'], label=batch, color=color, alpha=0.7)

    # Customize the plot
    
    plt.title('PCA Plot')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    #tsne plot
    
# tsne plot 
def perform_tsne_and_plot(df, random_state=seed):
    """
    Performs t-SNE on the dataset and plots the first two components.
    
    Parameters:
    df (pd.DataFrame): Input DataFrame with a 'batch' column.
    random_state (int): Random seed for reproducibility.
    
    Returns:
    None
    """
    # Extract features (exclude the 'batch' column)
    features = df.drop(columns=['Batch'])

    # Perform t-SNE
    tsne = TSNE(n_components=2, random_state=seed)
    tsne_components = tsne.fit_transform(features)

    # Create a DataFrame for the t-SNE results
    tsne_df = pd.DataFrame(data=tsne_components, columns=['TSNE1', 'TSNE2'])
    tsne_df['Batch'] = df['Batch'].values  # Add the batch column

    # Plot t-SNE results
    plt.figure(figsize=(8, 4))
    for batch, color in zip(['discovery', 'replication'], ['steelblue', 'coral']):
        subset = tsne_df[tsne_df['Batch'] == batch]
        plt.scatter(subset['TSNE1'], subset['TSNE2'], label=batch, color=color, alpha=0.7)

    # Customize the plot
    plt.title('t-SNE of Discovery and Replication Data')
    plt.xlabel('t-SNE Component 1')
    plt.ylabel('t-SNE Component 2')
    plt.legend()
    plt.grid(linestyle='--', linewidth=0.5, alpha=0.7)
    plt.show()


