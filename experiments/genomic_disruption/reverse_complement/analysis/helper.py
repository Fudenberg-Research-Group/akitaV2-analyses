import pandas as pd
from akita_utils.format_io import h5_to_df


def read_reverse_complement_disruption(disruption_dir, rc_disruption_dir):
    """
    Reads and processes disruption data for multiple models and their reverse complements.

    This function reads disruption data for models and their reverse complements from the provided directories.
    It calculates and returns the average values of specified metrics across all models.

    Parameters:
    disruption_dir (str): The directory path containing the disruption H5 files for the models.
    rc_disruption_dir (str): The directory path containing the reverse complement disruption H5 files for the models.

    Returns:
    pandas.DataFrame: A DataFrame with the disruption data, including average values for each metric across all models.
    """
    models = ["model_0", "model_1", "model_2", "model_3"]
    metrics = ["SCD", "INS-16", "INS-64"]
    
    # Initialize a dataframe to store all models' data
    df_all = pd.DataFrame()
    
    for i, model in enumerate(models):
        # Read model data
        df_model = h5_to_df(f"{disruption_dir}{model}.h5", metrics, average=True)
        df_rc_model = h5_to_df(f"{rc_disruption_dir}{model}_rc.h5", metrics, average=True)
        
        # Rename columns
        df_model = df_model.rename(columns={metric: f"{metric}_m{i}" for metric in metrics})
        df_rc_model = df_rc_model.rename(columns={metric: f"RC_{metric}_m{i}" for metric in metrics})
        
        # Concatenate data
        if df_all.empty:
            df_all = pd.concat([df_model, df_rc_model], axis=1)
        else:
            df_all = pd.concat([df_all, df_model, df_rc_model], axis=1)
    
    # Calculate averages
    for metric in metrics:
        df_all[metric] = df_all[[f"{metric}_m{i}" for i in range(4)]].mean(axis=1)
        df_all[f"RC_{metric}"] = df_all[[f"RC_{metric}_m{i}" for i in range(4)]].mean(axis=1)
    
    return df_all
