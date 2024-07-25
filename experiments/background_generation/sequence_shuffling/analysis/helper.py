import pandas as pd
from akita_utils.format_io import h5_to_df


def read_shuffling_data(data_dir, stat_names=["SCD"]):
    """
    Reads and processes shuffling data for multiple models.

    This function reads shuffling data from H5 files for a set of models, processes the specified statistics,
    and returns a DataFrame with the average values of these statistics across all models.

    Parameters:
    data_dir (str): The directory path containing the shuffling H5 files for the models.
    stat_names (list of str): The names of the statistics to be processed. Default is ["SCD"].

    Returns:
    pandas.DataFrame: A DataFrame with the shuffling data, including average values for each specified statistic across all models.
    """
    models = ["model_0", "model_1", "model_2", "model_3"]
    df_all = pd.DataFrame()

    for i, model in enumerate(models):
        df_model = h5_to_df(f"{data_dir}{model}.h5", stat_names, average=True)
        df_model = df_model.rename(columns={"SCD": f"SCD_m{i}"})

        if df_all.empty:
            df_all = df_model
        else:
            df_all = pd.concat([df_all, df_model[f"SCD_m{i}"]], axis=1)

    df_all["SCD"] = df_all[[f"SCD_m{i}" for i in range(4)]].mean(axis=1)

    return df_all
