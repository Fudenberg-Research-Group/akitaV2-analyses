import pandas as pd
from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import average_stat_over_targets


def read_multi_model_single_mutagenesis_data(
    data_dir, keys_to_ignore=["disruption_SCD", "insertion_SCD"]
):
    """
    Reads and processes mutagenesis data from multiple models stored in HDF5 files.

    For each model (0 to 3), this function:
    1. Loads data from an HDF5 file.
    2. Computes statistics for specified mutations (SCD, INS-16, INS-64).
    3. Aggregates these statistics into a single DataFrame.

    The final DataFrame includes the average statistics across all models for each mutation type.
    The function also drops unnecessary columns that include specific patterns in their names.

    Parameters:
    - data_dir (str): The directory containing the HDF5 files for each model.
    - keys_to_ignore (list of str): List of keys to ignore when loading the HDF5 data. Default is ["disruption_SCD", "insertion_SCD"].

    Returns:
    - pd.DataFrame: A DataFrame with the average statistics for each mutation type across all models, and with unnecessary columns removed.
    """

    stats = ["SCD", "INS-16", "INS-64"]
    df = pd.DataFrame()  # Initialize an empty DataFrame

    for model_index in range(4):
        model_df = h5_to_df(
            f"{data_dir}/model_{model_index}.h5",
            stats,
            average=False,
            ignore_keys=keys_to_ignore,
        )
        if model_index == 0:
            df = model_df[
                [
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "disruption_SCD",
                    "insertion_SCD",
                    "mutated_nucleotide",
                    "original_nucleotide",
                    "position",
                ]
            ].copy()
        for stat in stats:
            df_stat = average_stat_over_targets(
                model_df, model_index=model_index, head_index=1, stat=stat
            )
            df[f"{stat}_m{model_index}"] = df_stat[f"{stat}_m{model_index}"]

    for stat in stats:
        df[stat] = df[
            [f"{stat}_m{model_index}" for model_index in range(4)]
        ].mean(axis=1)

    return df
