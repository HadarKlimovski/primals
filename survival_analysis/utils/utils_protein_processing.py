import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler


from sklearn.preprocessing import StandardScaler

from lifelines import statistics

from configs.constants import  KERATIN,IMMUNE_PROTEINS,PROTEOMIC_DATA_REPLICATION,CLINICAL_DATA_DISCOVERY,CLINICAL_DATA_REP,DELTA_ALSFRS_REP,DELTA_ALSFRS_DIS,MAIN_DIR



# ----------------------------------------------------------------------------------------
def handle_duplicate_columns(df):
    """
    Standardizes column names by removing numeric suffixes and averages columns with identical names.

    Parameters:
    df (pd.DataFrame): The input DataFrame with potential duplicate columns.

    Returns:
    pd.DataFrame: A DataFrame with duplicate columns averaged.
    """
    # Remove numeric suffixes like '.1', '.2', etc., from column names
    df.columns = df.columns.str.replace(r'\.\d+$', '', regex=True)
    
    # Group by column names and calculate the mean for duplicates
    df = df.groupby(df.columns, axis=1).mean()
    
    return df

     

#plots thast describe the data
def describe_proteins_and_samples(df):
   
    # Use describe() to calculate descriptive statistics for proteins (columns)
    # protein_stats = df.describe()
    # print(protein_stats)
    
    # Extract mean and variance for histogram plotting
    protein_means = df.mean()
    protein_variances = df.var()
    

    # Plot the distribution of protein means
    plt.figure(figsize=(8, 4))
    plt.hist(protein_means, bins=20, edgecolor='black', alpha=0.7)
    plt.title('Distribution of Protein Means')
    plt.xlabel('Mean Protein Abundance')
    plt.ylabel('Frequency')
    plt.show()

    # Plot the distribution of protein variances
    plt.figure(figsize=(8, 4))
    plt.hist(protein_variances, bins=20, edgecolor='black', alpha=0.7)
    plt.title('Distribution of Protein Variances')
    plt.xlabel('Variance of Protein Abundance')
    plt.ylabel('Frequency')
    plt.show()

    # Calculate statistics for samples (rows)
    sample_protein_counts = df.count(axis=1)

    # Plot the distribution of protein counts per sample
    plt.figure(figsize=(8, 4))
    plt.hist(sample_protein_counts, bins=20, edgecolor='black', alpha=0.7)
    plt.title('Distribution of Protein Counts per Sample')
    plt.xlabel('Number of Detected Proteins per Sample')
    plt.ylabel('Frequency of Samples')
    plt.show()
    
    
    return df

def describe_clinical_data(df):
     # Plotting the box plot for "Survival_from_onset (months)"
    plt.figure(figsize=(8, 4))
    plt.boxplot(df["Survival_from_onset (months)"])
    plt.title("Box Plot of Survival from Onset (months)")
    plt.ylabel("Survival from Onset (months)")
    plt.grid(True)
    
    plt.figure(figsize=(8, 4))
    plt.hist(df["Survival_from_onset (months)"], bins=10, edgecolor='black', alpha=0.7)
    plt.title('Distribution of Survival time from Onset')
    plt.xlabel('Survival time(month)')
    plt.ylabel('Frequency of Samples')
    plt.show()
    
    plt.figure(figsize=(8, 4))
    df.boxplot(column="Survival_from_onset (months)", by="Sex", grid=False)
    num_males = df[df['Sex'] == 0].shape[0]
    num_females = df[df['Sex'] == 1].shape[0]
    plt.text(1, df["Survival_from_onset (months)"].max(), f'N = {num_males}', horizontalalignment='center', verticalalignment='center')
    plt.text(2, df["Survival_from_onset (months)"].max(), f'N = {num_females}', horizontalalignment='center', verticalalignment='center')
    
    plt.title("Box Plot of Survival from Onset (months) by Sex")
    plt.suptitle("")
    plt.xlabel("Sex")
    plt.ylabel("Survival from Onset (months)")
    plt.xticks([1, 2], ["Male", "Female"])  # Set custom x-tick labels
    plt.grid(True)
    plt.show()
    
    # Histogram for "Age Onset (years)" with average line
    plt.figure(figsize=(8, 4))
    plt.hist(df["Age Onset (years)"], bins=10, edgecolor='black', alpha=0.7)
    average_age_onset = df["Age Onset (years)"].mean()
    plt.axvline(average_age_onset, color='red', linestyle='dotted', linewidth=2)
    #plt.text(average_age_onset, plt.ylim()[1] * 0.9, f'Average: {average_age_onset:.2f} years', color='red', horizontalalignment='right')
    plt.title('Distribution of Age Onset (years)')
    plt.xlabel('Age Onset (years)')
    plt.ylabel('Frequency of Samples')
    plt.show()
    
    return df
# ----------------------------------------------------------------------------------------
#function to drop IMMUNE_PROTEINS and KERATIN proteins and proteins with low variance 

def drop_protein(df):
    """
    Parameters:
    IMMUNE_PROTEINS (list): List of immune protein columns to drop.
    KERATIN (list): List of keratin columns to drop.
    threshold (float): Variance threshold for filtering columns.
    """
    # Drop immune protein columns and print the shape
    df = df.drop(columns=IMMUNE_PROTEINS, errors='ignore')
    print(f"Shape after dropping immune proteins: {df.shape}")

    # Drop keratin columns and print the shape
    df = df.drop(columns=KERATIN, errors='ignore')
    print(f"Shape after dropping keratin: {df.shape}")

    # Remove columns with unrecognized names (e.g., NaN or None)
    valid_columns = df.columns.dropna()  # Drop columns where the name is NaN or None
    df = df[valid_columns]
    print(f"Shape after removing unrecognized column names: {df.shape}")

    # Filter out proteins that are found in less than 50% of the tested samples
    required_count = len(df) * 0.5
    df = df.loc[:, df.count() >= required_count]
    print(f"Shape after filtering proteins >50%: {df.shape}")

    return df





def low_variance_proteins(df, threshold=0.5):
    """
    Removes proteins (columns) with variance below the specified threshold.

    Parameters:
    df (pd.DataFrame): DataFrame containing protein abundance data.
    threshold (float): Variance threshold for filtering proteins.

    Returns:
    pd.DataFrame: Filtered DataFrame with proteins having variance above the threshold.
    """
    # Calculate variance and apply variance threshold filtering
    selector = VarianceThreshold(threshold=threshold)
    selector.fit(df)

    # Get the mask for proteins with variance above the threshold
    mask = selector.get_support()
    
    # Filter the DataFrame using the mask to avoid duplication
    df_filtered = df.loc[:, mask]
    print(f"Proteins with variance above threshold (threshold={threshold}): {df_filtered.shape}")

    # Plot the distribution of variances with a threshold line
    variances = df.var()  # Variance of original DataFrame
    plt.figure(figsize=(8, 4))
    plt.hist(variances, bins=40, alpha=0.7, edgecolor='black')
    plt.axvline(x=threshold, color='red', linestyle='--', label=f'Threshold = {threshold}')
    plt.title('Distribution of Protein Variances with Threshold Line')
    plt.xlabel('Variance of Protein Abundance')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)
    plt.show()

    return df_filtered
#-----


def identify_low_variance_and_high_correlation(df, variance_threshold=0.05, correlation_threshold=0.8):
    """
    Identifies and prints proteins with low variance and pairs of proteins with high correlation.

    Parameters:
    df (DataFrame): DataFrame containing protein abundance data.
    variance_threshold (float): Threshold below which a protein is considered to have low variance.
    correlation_threshold (float): Threshold above which proteins are considered highly correlated.

    Returns:
    Tuple containing:
    - List of low variance proteins.
    - DataFrame of highly correlated protein pairs with their correlation values.
    """
    # Step 1: Identify low variance proteins
    selector = VarianceThreshold(threshold=variance_threshold)
    selector.fit(df)
    low_variance_proteins = df.columns[~selector.get_support()].tolist()

    print(f"Proteins with variance below {variance_threshold}:")
    print(low_variance_proteins)

    # Step 2: Identify highly correlated proteins
    correlation_matrix = df.corr(method='spearman').abs()
    upper_triangle = np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool)
    upper_triangle_corr = correlation_matrix.where(upper_triangle)

    correlated_pairs = []
    for column in upper_triangle_corr.columns:
        high_corr = upper_triangle_corr[column][upper_triangle_corr[column] > correlation_threshold]
        for idx, corr_value in high_corr.items():
            correlated_pairs.append((column, idx, corr_value))

    correlated_df = pd.DataFrame(correlated_pairs, columns=['Protein_1', 'Protein_2', 'Correlation'])
    print(f"Pairs of proteins with Spearman correlation above {correlation_threshold}:")
    print(correlated_df)

    return low_variance_proteins, correlated_df

import pandas as pd

def drop_unnecessary_columns(df, columns_to_drop):
    """
    Removes specified columns from a DataFrame.

    Parameters:
    df (pd.DataFrame): The DataFrame from which columns will be removed.
    columns_to_drop (list): A list of column names to be dropped.

    Returns:
    pd.DataFrame: A DataFrame with the specified columns removed.
    """
    # Drop the specified columns
    df = df.drop(columns=columns_to_drop, errors='ignore')
    return df



# --------------------------------------------------------------------------------

def norm_log(df, apply_norm=False):

    # Apply log2 transformation, adding a small value to avoid log(0)
    df = df.astype(float)
    df = np.log2(df + 1e-10)  

    if apply_norm:
         # Normalize using L1 normalization
        row_sums = df.abs().sum(axis=1)
        df = df.div(row_sums, axis=0)
         # Check that the sum of each row is 1
        row_sums_check = df.sum(axis=1)
        print(f"Sum of rows is 1: {np.allclose(row_sums_check, 1)}")
        print(f"Sum of the first 3 rows after normalization: {row_sums_check.head(3)}")

    return df

# ----------------------------------------------------------------------------
#Merg protein data with clinical data.

def merge_proteins_and_clinical_data_discovery(df):
    clinical_data = pd.read_csv(CLINICAL_DATA_DISCOVERY)
    print(f"Loading clinical data from {CLINICAL_DATA_DISCOVERY}\nShape: {clinical_data.shape}")
    
    clinical_data1 = clinical_data[["sample number", "ALSFRS score (unit)","Age Onset (years)",'Sex','Disease Format']]
    clinical_data1.set_index("sample number", inplace=True)
    df_with_clinical = df.join(clinical_data1, how='inner')
    # scaler = StandardScaler()
    # df_with_clinical = pd.DataFrame(scaler.fit_transform(df_with_clinical), 
    #                                columns=df_with_clinical.columns, 
    #                                index=df_with_clinical.index)
    clinical_data2 = clinical_data[["sample number","Survival_from_onset (months)","Status dead=1"]]
    clinical_data2.set_index("sample number", inplace=True)
    df_with_clinical = df_with_clinical.join(clinical_data2, how='inner')
    
    sex_mapping = {'M': 0, 'F': 1}
    df_with_clinical['Sex'] = df_with_clinical['Sex'].map(sex_mapping)
    disease_format_mapping = {'Limb': 0, 'Bulbar': 1}
    df_with_clinical['Disease Format'] = df_with_clinical['Disease Format'].map(disease_format_mapping)
    return df_with_clinical


def merge_proteins_and_clinical_data_rep(df):
    clinical_data = pd.read_csv(CLINICAL_DATA_REP)
    print(f"Loading clinical data from {CLINICAL_DATA_REP}\nShape: {clinical_data.shape}")
    
    clinical_data1 = clinical_data[["sample number", "ALSFRS score (unit)","Age Onset (years)",'Sex','Disease Format']]
    clinical_data1.set_index("sample number", inplace=True)
    df_with_clinical = df.join(clinical_data1, how='inner')
    # scaler = StandardScaler()
    # df_with_clinical = pd.DataFrame(scaler.fit_transform(df_with_clinical), 
    #                                columns=df_with_clinical.columns, 
    #                                index=df_with_clinical.index)
    clinical_data2 = clinical_data[["sample number","Survival_from_onset (months)","Status dead=1"]]
    clinical_data2.set_index("sample number", inplace=True)
    df_with_clinical = df_with_clinical.join(clinical_data2, how='inner')
    
    sex_mapping = {'M': 0, 'F': 1}
    df_with_clinical['Sex'] = df_with_clinical['Sex'].map(sex_mapping)
    disease_format_mapping = {'Limb': 0, 'Bulbar': 1}
    df_with_clinical['Disease Format'] = df_with_clinical['Disease Format'].map(disease_format_mapping)
    return df_with_clinical



def scale_data(df, scaling="zscore"):
    
    # Identify protein columns (excluding clinical data)
    protein_columns = [col for col in df.columns if col not in ['Survival_from_onset (months)', 'Status dead=1']]
    
    # Apply the requested scaling method
    if scaling == "zscore":
        scaler = StandardScaler()
        df[protein_columns] = scaler.fit_transform(df[protein_columns])
        print("Z-score normalization applied.")
    elif scaling == "minmax":
        scaler = MinMaxScaler()
        df[protein_columns] = scaler.fit_transform(df[protein_columns])
        print("Min-Max scaling applied.")
    else:
        print("No scaling applied.")
    
    return df


def plot_combined_distributions(data, df_imputed):
    # Calculate row-wise means for each sample
    df_means = data.mean(axis=1)  # Mean before imputation
    df_imputed_means = df_imputed.mean(axis=1)  # Mean after imputation

    # Create a combined DataFrame for plotting
    combined_data = pd.DataFrame({
        'Value': pd.concat([df_means, df_imputed_means], axis=0),
        'Data': ['Before Imputation'] * len(df_means) + ['After Imputation'] * len(df_imputed_means)
    })

    # Plot overlapping distributions using seaborn
    plt.figure(figsize=(8, 6))
    sns.kdeplot(data=combined_data, x='Value', hue='Data', fill=True, common_norm=False, alpha=0.5)
    plt.title('Distribution of Imputed and Non-Imputed Data')
    plt.xlabel('Mean Protein Abundance')
    plt.ylabel('Density')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()



def plot_nan_histogram(df):
    # Calculate the number of NaN values per sample (row)
    nan_counts = df.isna().sum(axis=1)
    
    # Plot the histogram
    plt.figure(figsize=(8, 6))
    plt.hist(nan_counts, bins=20, edgecolor='black', alpha=0.7)
    plt.title('Histogram of NaN Counts per Sample')
    plt.xlabel('Number of NaN Values')
    plt.ylabel('Frequency (Number of Samples)')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
