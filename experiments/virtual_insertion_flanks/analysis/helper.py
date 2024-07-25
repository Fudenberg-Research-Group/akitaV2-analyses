import pandas as pd
from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import average_stat_over_targets


def read_multi_model_single_flanks_data(data_dir):
    """
    Reads and processes data from multiple models and backgrounds, then computes average statistics.

    This function aggregates data from multiple models and background conditions to compute average statistics
    for different features. It reads HDF5 files for each model and background, processes the data to average
    statistics over targets, and combines the results into a single DataFrame.

    The function performs the following steps:
    1. Reads and processes data for model 0 and background 0 to initialize the DataFrame.
    2. Iteratively reads and processes data for other models and backgrounds.
    3. Computes and includes average statistics for "SCD", "INS-16", and "INS-64" across all models and backgrounds.
    4. Computes overall averages for these statistics across all models and backgrounds.

    Args:
        data_dir (str): The directory containing the HDF5 files for the models and backgrounds. The directory
                        should have subdirectories for each model, and within each model subdirectory, 
                        there should be files for each background.

    Returns:
        pd.DataFrame: A DataFrame containing the averaged statistics for "SCD", "INS-16", and "INS-64" across
                      all models and backgrounds, along with other relevant columns.
    """
    def process_data(model_index, bg_index, stat):
        """Helper function to process and rename data."""
        df = h5_to_df(f"{data_dir}/model_{model_index}/model_{model_index}_bg_{bg_index}.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
        df_avg = average_stat_over_targets(df, model_index=model_index, head_index=1, stat=stat)
        return df_avg.rename(columns={f"{stat}_m{model_index}": f"{stat}_m{model_index}_bg{bg_index}"})[f"{stat}_m{model_index}_bg{bg_index}"]

    # Read and process data for model 0, background 0 only
    df_m0_bg0 = h5_to_df(f"{data_dir}/model_0/model_0_bg_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    # df = pd.DataFrame()
    df = df_m0_bg0[["chrom", "start", "end", "strand", "flank_bp", "insertion_SCD", "disruption_SCD"]].copy()
    for stat in ["SCD", "INS-16", "INS-64"]:
        df_stat = average_stat_over_targets(df_m0_bg0, model_index=0, head_index=1, stat=stat)
        df[f"{stat}_m0_bg0"] = df_stat[f"{stat}_m0"]

    # Process data for other models and backgrounds
    for model_index in range(4):
        start_from = 1 if model_index == 0 else 0
        for bg_index in range(start_from, 10):
            for stat in ["SCD", "INS-16", "INS-64"]:
                df[f"{stat}_m{model_index}_bg{bg_index}"] = process_data(model_index, bg_index, stat)

    # Averaging
    for stat in ["SCD", "INS-16", "INS-64"]:
        df[stat] = df[[f"{stat}_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)

    return df


def read_multi_model_double_flanks_data(data_dir):
    """
    Reads and processes data from multiple models with double flank conditions, then computes average statistics.

    This function aggregates data from multiple models with double flank conditions to compute average statistics
    for different features. It reads HDF5 files for each model with double flanks, processes the data to average
    statistics over targets, combines the results, and computes overall averages.

    The function performs the following steps:
    1. Reads and concatenates data from two files for each model (with double flanks).
    2. Processes the concatenated data to compute average statistics for "SCD", "INS-16", and "INS-64".
    3. Drops columns related to individual targets.
    4. Combines data from all models into a single DataFrame.
    5. Aggregates the data by grouping based on specific columns and computing mean statistics for each group.
    6. Computes overall averages for the aggregated statistics across all models.

    Args:
        data_dir (str): The directory containing the HDF5 files for the models with double flank conditions. 
                        The directory should contain files for each model, with filenames following the pattern 
                        "model_<model_index>_<flank_index>.h5".

    Returns:
        pd.DataFrame: A DataFrame containing the averaged statistics for "SCD", "INS-16", and "INS-64" across
                      all models, along with other relevant columns. The DataFrame includes mean values computed
                      for each group defined by "chrom", "start", "end", "strand", "orientation", and "flank_bp".
    """
    df_combined = None

    for model_index in range(4):
        df_0 = h5_to_df(f"{data_dir}model_{model_index}_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
        df_1 = h5_to_df(f"{data_dir}model_{model_index}_1.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
        df = pd.concat([df_0, df_1])
        df = average_stat_over_targets(df, model_index=model_index, head_index=1, stat="SCD")
        df[f"INS-16_m{model_index}"] = average_stat_over_targets(df, model_index=model_index, head_index=1, stat="INS-16")[f"INS-16_m{model_index}"]
        df[f"INS-64_m{model_index}"] = average_stat_over_targets(df, model_index=model_index, head_index=1, stat="INS-64")[f"INS-64_m{model_index}"]
        df = df.drop(columns=[f"SCD_h1_m{model_index}_t{target_ind}" for target_ind in range(6)] +
                         [f"INS-16_h1_m{model_index}_t{target_ind}" for target_ind in range(6)] +
                         [f"INS-64_h1_m{model_index}_t{target_ind}" for target_ind in range(6)])

        if df_combined is None:
            df_combined = df
        else:
            df_combined[f"SCD_m{model_index}"] = df[f"SCD_m{model_index}"]
            df_combined[f"INS-16_m{model_index}"] = df[f"INS-16_m{model_index}"]
            df_combined[f"INS-64_m{model_index}"] = df[f"INS-64_m{model_index}"]

    df = df_combined.groupby(["chrom", "start", "end", "strand", "orientation", "flank_bp"]).agg({
        "SCD_m0": "mean",
        "SCD_m1": "mean",
        "SCD_m2": "mean",
        "SCD_m3": "mean",
        "INS-16_m0": "mean",
        "INS-16_m1": "mean",
        "INS-16_m2": "mean",
        "INS-16_m3": "mean",
        "INS-64_m0": "mean",
        "INS-64_m1": "mean",
        "INS-64_m2": "mean",
        "INS-64_m3": "mean"
    }).reset_index()

    df["SCD"] = df[[f"SCD_m{model_ind}" for model_ind in range(4)]].mean(axis=1)
    df["INS-16"] = df[[f"INS-16_m{model_ind}" for model_ind in range(4)]].mean(axis=1)
    df["INS-64"] = df[[f"INS-64_m{model_ind}" for model_ind in range(4)]].mean(axis=1)

    return df