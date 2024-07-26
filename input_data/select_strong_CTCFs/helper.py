import pandas as pd
from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import (
    average_stat_over_targets,
    average_stat_over_backgrounds,
)


def read_and_average_virtual_exp(
    data_dir,
    stat_to_average="SCD",
    all_calculated_stats=["SCD", "INS-16", "INS-64"],
    model_numbers=8,
    head_index=1,
):
    """
    Reads data from h5 files for different models, calculates the average of a specified statistic over various targets and backgrounds, and then averages these values over all models.

    Parameters:
    - data_dir (str): The directory path where the h5 files are stored.
    - stat_to_average (str, optional): The specific statistic to be averaged. Defaults to "SCD".
    - all_calculated_stats (list of str, optional): List of all statistics calculated in the data. Defaults to ["SCD", "INS-16", "INS-64"].
    - model_numbers (int, optional): The number of models to consider for averaging. Defaults to 8.

    Returns:
    - DataFrame: A pandas DataFrame containing the averaged statistics for each chromosomal location across all models.
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
            df, model_index=i, head_index=head_index, stat=stat_to_average
        )
        for i, df in enumerate(dfs)
    ]

    print("averaging over backgrounds")
    df_tgbg = [
        average_stat_over_backgrounds(
            df, model_index=i, head_index=head_index, stat=stat_to_average
        )
        for i, df in enumerate(df_tg)
    ]

    print(f"collecting data for {stat_to_average}")

    # Start with the first DataFrame and add the statistics columns
    stat_collected = df_tgbg[0][["chrom", "end", "start", "strand"]].copy()

    for i in range(model_numbers):
        stat_collected[f"{stat_to_average}_m{i}"] = df_tgbg[i][
            f"{stat_to_average}_m{i}"
        ]

    # Ensure all columns are numeric and fill NaNs if needed
    stat_collected = stat_collected.fillna(
        0
    )  # You can use other strategies like forward/backward fill if needed

    # Averaging over models
    stat_collected[stat_to_average] = stat_collected[
        [f"{stat_to_average}_m{i}" for i in range(model_numbers)]
    ].mean(axis=1)

    return stat_collected


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
