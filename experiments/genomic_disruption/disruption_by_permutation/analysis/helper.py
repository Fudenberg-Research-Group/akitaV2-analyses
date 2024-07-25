import pandas as pd
from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import average_stat_over_targets


def read_and_average_genomic_exp(
    data_dir,
    stat_to_average="SCD",
    all_calculated_stats=["SCD", "INS-16", "INS-64"],
    model_numbers=8,
):
    """
    Reads data from h5 files and calculates the average of a specified statistic across multiple models.

    Parameters:
    - data_dir (str): The directory path where the h5 files are stored.
    - stat_to_average (str, optional): The statistic to be averaged. Default is "SCD".
    - all_calculated_stats (list of str, optional): A list of all statistics calculated in the data. Default is ["SCD", "INS-16", "INS-64"].
    - model_numbers (int, optional): The number of models to consider for averaging. Default is 8.

    Returns:
    - DataFrame: A pandas DataFrame containing the concatenated data from all models with the average of the specified statistic.
    """

    print("reading h5 files to dataframes")
    dfs = [
        h5_to_df(
            f"{data_dir}/model_{i}.h5", all_calculated_stats, average=False
        )
        for i in range(model_numbers)
    ]

    print("averaging over targets")
    df_tg = [
        average_stat_over_targets(
            df, model_index=i, head_index=1, stat=stat_to_average
        )
        for i, df in enumerate(dfs)
    ]

    print(f"collecting data for {stat_to_average}")
    # Collect data from all models
    data_frames = [
        df_tg[i][["chrom", "end", "start", "strand"]].assign(
            **{f"{stat_to_average}_m{i}": df_tg[i][f"{stat_to_average}_m{i}"]}
        )
        for i in range(model_numbers)
    ]

    # Concatenate all DataFrames
    stat_collected = pd.concat(data_frames, axis=1)

    # Averaging over models
    stat_collected[stat_to_average] = stat_collected[
        [f"{stat_to_average}_m{i}" for i in range(model_numbers)]
    ].mean(axis=1)

    return stat_collected


def read_disruption_smf_data(data_dir):
    """
    Reads and processes disruption data for multiple models, averaging specified statistics.

    Parameters:
    - data_dir (str): The directory path where the model H5 files are stored.

    Returns:
    - DataFrame: A pandas DataFrame containing the averaged statistics for all models.
    """
    # List of models
    model_indices = range(4)

    # List to store dataframes for each model
    dfs = []

    # Read and process data for each model
    for i in model_indices:
        df = h5_to_df(
            f"{data_dir}model_{i}.h5",
            ["SCD", "INS-16", "INS-64"],
            average=True,
        )
        df = df.rename(
            columns={
                "SCD": f"SCD_m{i}",
                "INS-16": f"INS-16_m{i}",
                "INS-64": f"INS-64_m{i}",
            }
        )
        if i == 0:
            df = df.drop(
                columns=["isBound", "name", "TF", "chipseq.score", "width"]
            )
        dfs.append(df)

    # Merge all model DataFrames on the rownames column, with appropriate suffixes
    df_merged = dfs[0]
    for i in range(1, len(dfs)):
        df_merged = df_merged.merge(
            dfs[i], on="rownames", suffixes=("", f"_m{i}")
        )

    # Rename the rownames column to TFBS_cluster
    df_merged = df_merged.rename(columns={"rownames": "TFBS_cluster"})

    # Calculate the average for each statistic across models
    for stat in ["SCD", "INS-16", "INS-64"]:
        df_merged[stat] = df_merged[
            [f"{stat}_m{i}" for i in model_indices]
        ].mean(axis=1)

    # Drop duplicates based on the TFBS_cluster column
    df_merged.drop_duplicates("TFBS_cluster", inplace=True)

    return df_merged
