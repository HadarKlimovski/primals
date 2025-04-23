from lifelines import KaplanMeierFitter, NelsonAalenFitter, CoxPHFitter
import matplotlib.pyplot as plt
 

def cox_proportional_hazard_model(data,
                                  survival_duration,
                                  survival_status,
                                  model,
                                  strata=None,
                                  covariate_groups=None,
                                  show_plot=True):
    """
    CPH regression model
       
        # alpha - confidence intervals
 
        # penalizer - an L2 penalizer to the size of the coefficients during regression.
            This improves stability of the estimates and controls for high correlation between covariates.
            For example, this shrinks the absolute value of 𝛽𝑖. The penalty is 12penalizer||𝛽||2.
       
        # strata - list of columns to use in stratification. Useful if a categorical covariate does not obey the proportional hazard assumption.
    ### higher cv, lower li_ratio, and lower penalizer shrinks the pval = more significant results
    """
   
    cph = CoxPHFitter(alpha=0.05, penalizer=0.01, strata=None, l1_ratio=0.0)
   
    cph = cph.fit(
                df=data,
                duration_col=survival_duration,
                event_col=survival_status,
                strata=strata,
                show_progress=False)
    ##cph.check_assumptions(data, p_value_threshold=0.05, show_plots=True)
 
   
    print(
        "\n\n Concordance index of the model", cph.concordance_index_)#,
        #"\n baseline hazard", cph.baseline_hazard_,
        #"\n baseline cumulative hazard_", cph.baseline_cumulative_hazard_,
        #"\n baseline survival", cph.baseline_survival_,
        #"\n the variance matrix of the coefficients", cph.variance_matrix_ ,
        #"\n\n The estimated coefficients \n", cph.hazards_)
 
    
    if data.shape[1] <= 11:
        fig_size = (10, 6)
        axis_label_size = 12
        title_size = 18
    else:
        fig_size = (20, 12)
        axis_label_size = 10
        title_size = 18
       
    # ---------------------------------------------------------------------
    # forest plot
    if show_plot:
        fig, ax = plt.subplots(figsize=(fig_size[0], fig_size[1]))
        cph.plot(c='b', ax=ax)
    
        # set axis label names and title
        #ax.tick_params(labelsize=axis_label_size, rotation=20, labelcolor='#020381')
        # set your ticks manually
        plt.title("Cox Proportional Hazards Regression Analysis - Concordance-index: " + str(round(cph.concordance_index_,2)), y=1.05, size=26, color='#960203', family='fantasy')
    
        # present the hazard value on the plot
        for i, v in enumerate(cph.params_.sort_values()):
            ax.text(v - 0.05, i + .1, str(round(v, 2)), color='#960203', alpha=0.9, size=axis_label_size, fontweight='bold', rotation=10, family='fantasy')
    
    # ---------------------------------------------------------------------
    if show_plot:
        cph.print_summary(model=model)
    
    # ---------------------------------------------------------------------
    if not covariate_groups is None:
        cph.plot_covariate_groups(covariates=covariate_groups, values=[0,1])
   
    #-----------------------------------------------------------------
    ## Lets predict the survival curve for the selected customers.
    ## Customers can be identified with the help of the number mentioned against each curve.
    ##cph.predict_survival_function(data).plot()
 
    return cph


def dichotomize_proteins_by_median(df, protein_columns):
    """
    Convert selected protein columns to dichotomous (binary) format based on their median values.

    :param df: DataFrame containing protein expression data.
    :param protein_columns: List of protein column names to dichotomize.
    :return: A new DataFrame containing only the dichotomized binary columns along with the specified unchanged columns.
    """
    # Ensure required columns exist
    missing_cols = [col for col in protein_columns if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Columns not found in DataFrame: {missing_cols}")

    df_binary = df.copy()

    # Dichotomize the protein columns
    for protein in protein_columns:
        median_value = df_binary[protein].median()
        df_binary[protein] = (df_binary[protein] > median_value).astype(int)  # Replace with binary values

    # Ensure these columns remain unchanged at the end
    unchanged_columns = ['Survival_from_onset (months)', 'Status dead=1']
    final_columns = [col for col in protein_columns if col in df_binary.columns] + unchanged_columns

    # Return only the dichotomized columns and the unchanged ones at the end
    return df_binary[final_columns]

