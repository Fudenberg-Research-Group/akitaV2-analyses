import numpy as np
import pandas as pd
import pysam
import scipy

from akita_utils.format_io import h5_to_df


########################################
# Averaging functions                  #
########################################


def average_stat_over_targets(df, model_index, head_index, stat="SCD"):
    """
    Calculate the average of a specified statistical metric (stat) over multiple targets for a given model and head.

    Parameters:
    df (DataFrame): The input DataFrame containing the data.
    model_index (int): The index of the model for which the metric is calculated.
    head_index (int): The index of the head for which the metric is calculated.
    stat (str, optional): The statistical metric to calculate the average for (default is "SCD").

    Returns:
    DataFrame: A DataFrame with a new column containing the average of the specified metric for the specified model and head.
    """
    if head_index == 1:
        target_indices = 6
    else:
        target_indices = 5

    df[f"{stat}_m{model_index}"] = df[
        [
            f"{stat}_h{head_index}_m{model_index}_t{target_index}"
            for target_index in range(target_indices)
        ]
    ].mean(axis=1)
    return df


# def organize_df_by_target(input_df, stat, num_targets=6, head_index=1, model_index=0):
#     """
#     Extracts data from the specified columns in the input DataFrame to create a new DataFrame
#     with 'target' indices and corresponding values for the given statistical measure.

#     Parameters:
#     - input_df (pd.DataFrame): The input DataFrame containing the data.
#     - stat (str): The base name of the statistical measure columns.
#     - num_targets (int, optional): The number of target indices. Default is 6.
#     - head_index (int, optional): The head index value. Default is 1.
#     - model_index (int, optional): The model index value. Default is 0.

#     Returns:
#     - pd.DataFrame: A new DataFrame with 'target' indices and corresponding values for the given statistical measure.
#     """
#     target_indices = []
#     stat_score = []

#     for target_index in range(num_targets):
#         column_name = f"{stat}_h{head_index}_m{model_index}_t{target_index}"
#         target_indices.extend([target_index] * len(input_df[column_name]))
#         stat_score.extend(list(input_df[column_name]))

#     output_df = pd.DataFrame({"target": target_indices, stat: stat_score})
#     return output_df


def average_stat_over_backgrounds(
    df,
    model_index=0,
    head_index=1,
    num_backgrounds=10,
    stat="SCD",
    columns_to_keep=["chrom", "end", "start", "strand", "seq_id"],
    keep_background_columns=True,
):
    """
    Calculate the average of a specified statistical metric (stat) over multiple background samples for a given model and head.

    Parameters:
    df (DataFrame): The input DataFrame containing the data, including background information.
    model_index (int, optional): The index of the model for which the metric is calculated (default is 0).
    head_index (int, optional): The index of the head for which the metric is calculated (default is 1).
    num_backgrounds (int, optional): The number of background samples to consider (default is 10).
    stat (str, optional): The statistical metric to calculate the average for (default is "SCD").
    columns_to_keep (list, optional): A list of columns to keep in the output DataFrame (default is ["chrom", "end", "start", "strand", "seq_id"]).
    keep_background_columns (bool, optional): Whether to keep individual background columns in the output DataFrame (default is True).

    Returns:
    DataFrame: A DataFrame with the specified statistical metric's average for the specified model and head, along with optional columns.
    """

    if head_index == 1:
        target_indices = 6
    else:
        target_indices = 5

    num_sites = len(df) // num_backgrounds
    output_df = df[columns_to_keep][:num_sites]

    for bg_index in range(num_backgrounds):
        output_df[f"{stat}_bg{bg_index}"] = df[
            df["background_index"] == bg_index
        ][f"{stat}_m{model_index}"].values

    output_df[f"{stat}_m{model_index}"] = output_df[
        [f"{stat}_bg{bg_index}" for bg_index in range(num_backgrounds)]
    ].mean(axis=1)

    if keep_background_columns == False:
        output_df = output_df.drop(
            columns=[
                f"{stat}_bg{bg_index}" for bg_index in range(num_backgrounds)
            ]
        )

    return output_df


def average_stat_for_shift(df, shift, model_index, head_index, stat="SCD"):
    if head_index == 1:
        target_indices = 6
    else:
        target_indices = 5

    df[f"{stat}_{shift}"] = df[
        [
            f"{stat}_h{head_index}_m{model_index}_t{target_index}"
            for target_index in range(target_indices)
        ]
    ].mean(axis=1)
    return df


########################################
# HF5 reading functions                #
########################################


def read_and_average_virtual_exp(
    data_dir,
    stat_to_average="SCD",
    all_calculated_stats=["SCD", "INS-16", "INS-64"],
    model_numbers=8,
):
    """
    Reads data from h5 files for different models, calculates the average of a specified statistic over various targets and backgrounds, and then averages these values over all models.

    This function processes data from multiple models stored in h5 files. It averages a specified statistic first over different targets, then over different backgrounds for each model, and finally averages these values across all models.

    Parameters:
    - data_dir (str): The directory path where the h5 files are stored.
    - stat_to_average (str, optional): The specific statistic to be averaged. Defaults to "SCD".
    - all_calculated_stats (list of str, optional): List of all statistics calculated in the data. Defaults to ["SCD", "INS-16", "INS-64"].
    - model_numbers (int, optional): The number of models to consider for averaging. Defaults to 8.

    Returns:
    - DataFrame: A pandas DataFrame containing the averaged statistics for each chromosomal location across all models.
    """
    print("reading h5 files to dataframes")
    df_m0 = h5_to_df(
        data_dir + "/model_0.h5", all_calculated_stats, average=False
    )
    df_m1 = h5_to_df(
        data_dir + "/model_1.h5", all_calculated_stats, average=False
    )
    df_m2 = h5_to_df(
        data_dir + "/model_2.h5", all_calculated_stats, average=False
    )
    df_m3 = h5_to_df(
        data_dir + "/model_3.h5", all_calculated_stats, average=False
    )
    df_m4 = h5_to_df(
        data_dir + "/model_4.h5", all_calculated_stats, average=False
    )
    df_m5 = h5_to_df(
        data_dir + "/model_5.h5", all_calculated_stats, average=False
    )
    df_m6 = h5_to_df(
        data_dir + "/model_6.h5", all_calculated_stats, average=False
    )
    df_m7 = h5_to_df(
        data_dir + "/model_7.h5", all_calculated_stats, average=False
    )

    print("averaging over targets")
    df_m0_tg = average_stat_over_targets(
        df_m0, model_index=0, head_index=1, stat=stat_to_average
    )
    df_m1_tg = average_stat_over_targets(
        df_m1, model_index=1, head_index=1, stat=stat_to_average
    )
    df_m2_tg = average_stat_over_targets(
        df_m2, model_index=2, head_index=1, stat=stat_to_average
    )
    df_m3_tg = average_stat_over_targets(
        df_m3, model_index=3, head_index=1, stat=stat_to_average
    )
    df_m4_tg = average_stat_over_targets(
        df_m4, model_index=4, head_index=1, stat=stat_to_average
    )
    df_m5_tg = average_stat_over_targets(
        df_m5, model_index=5, head_index=1, stat=stat_to_average
    )
    df_m6_tg = average_stat_over_targets(
        df_m6, model_index=6, head_index=1, stat=stat_to_average
    )
    df_m7_tg = average_stat_over_targets(
        df_m7, model_index=7, head_index=1, stat=stat_to_average
    )

    print("averaging over backgrounds")
    df_m0_tgbg = average_stat_over_backgrounds(
        df_m0_tg, model_index=0, head_index=1, stat=stat_to_average
    )
    df_m1_tgbg = average_stat_over_backgrounds(
        df_m1_tg, model_index=1, head_index=1, stat=stat_to_average
    )
    df_m2_tgbg = average_stat_over_backgrounds(
        df_m2_tg, model_index=2, head_index=1, stat=stat_to_average
    )
    df_m3_tgbg = average_stat_over_backgrounds(
        df_m3_tg, model_index=3, head_index=1, stat=stat_to_average
    )
    df_m4_tgbg = average_stat_over_backgrounds(
        df_m4_tg, model_index=4, head_index=1, stat=stat_to_average
    )
    df_m5_tgbg = average_stat_over_backgrounds(
        df_m5_tg, model_index=5, head_index=1, stat=stat_to_average
    )
    df_m6_tgbg = average_stat_over_backgrounds(
        df_m6_tg, model_index=6, head_index=1, stat=stat_to_average
    )
    df_m7_tgbg = average_stat_over_backgrounds(
        df_m7_tg, model_index=7, head_index=1, stat=stat_to_average
    )

    print(f"collecting data for {stat_to_average}")
    stat_collected = pd.concat(
        [
            df_m0_tgbg["chrom"],
            df_m0_tgbg["end"],
            df_m0_tgbg["start"],
            df_m0_tgbg["strand"],
            df_m0_tgbg[f"{stat_to_average}_m0"],
            df_m1_tgbg[f"{stat_to_average}_m1"],
            df_m2_tgbg[f"{stat_to_average}_m2"],
            df_m3_tgbg[f"{stat_to_average}_m3"],
            df_m4_tgbg[f"{stat_to_average}_m4"],
            df_m5_tgbg[f"{stat_to_average}_m5"],
            df_m6_tgbg[f"{stat_to_average}_m6"],
            df_m7_tgbg[f"{stat_to_average}_m7"],
        ],
        axis=1,
    )

    # averaging over models
    stat_collected[f"{stat_to_average}"] = stat_collected[
        [
            f"{stat_to_average}_m{model_index}"
            for model_index in range(model_numbers)
        ]
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

    This function is designed to process genomic data stored in h5 format. It reads data from specified h5 files for different models, averages the selected statistic over targets for each model, and then computes the overall average of this statistic across all models.

    Parameters:
    - data_dir (str): The directory path where the h5 files are stored.
    - stat_to_average (str, optional): The statistic to be averaged. Default is "SCD".
    - all_calculated_stats (list of str, optional): A list of all statistics calculated in the data. Default is ["SCD", "INS-16", "INS-64"].
    - model_numbers (int, optional): The number of models to consider for averaging. Default is 8.

    Returns:
    - DataFrame: A pandas DataFrame containing the concatenated data from all models with the average of the specified statistic.
    """

    print("reading h5 files to dataframes")
    df_m0 = h5_to_df(
        data_dir + "/model_0.h5", all_calculated_stats, average=False
    )
    df_m1 = h5_to_df(
        data_dir + "/model_1.h5", all_calculated_stats, average=False
    )
    df_m2 = h5_to_df(
        data_dir + "/model_2.h5", all_calculated_stats, average=False
    )
    df_m3 = h5_to_df(
        data_dir + "/model_3.h5", all_calculated_stats, average=False
    )
    df_m4 = h5_to_df(
        data_dir + "/model_4.h5", all_calculated_stats, average=False
    )
    df_m5 = h5_to_df(
        data_dir + "/model_5.h5", all_calculated_stats, average=False
    )
    df_m6 = h5_to_df(
        data_dir + "/model_6.h5", all_calculated_stats, average=False
    )
    df_m7 = h5_to_df(
        data_dir + "/model_7.h5", all_calculated_stats, average=False
    )

    print("averaging over targets")
    df_m0_tg = average_stat_over_targets(
        df_m0, model_index=0, head_index=1, stat=stat_to_average
    )
    df_m1_tg = average_stat_over_targets(
        df_m1, model_index=1, head_index=1, stat=stat_to_average
    )
    df_m2_tg = average_stat_over_targets(
        df_m2, model_index=2, head_index=1, stat=stat_to_average
    )
    df_m3_tg = average_stat_over_targets(
        df_m3, model_index=3, head_index=1, stat=stat_to_average
    )
    df_m4_tg = average_stat_over_targets(
        df_m4, model_index=4, head_index=1, stat=stat_to_average
    )
    df_m5_tg = average_stat_over_targets(
        df_m5, model_index=5, head_index=1, stat=stat_to_average
    )
    df_m6_tg = average_stat_over_targets(
        df_m6, model_index=6, head_index=1, stat=stat_to_average
    )
    df_m7_tg = average_stat_over_targets(
        df_m7, model_index=7, head_index=1, stat=stat_to_average
    )

    print(f"collecting data for {stat_to_average}")
    stat_collected = pd.concat(
        [
            df_m1_tg["chrom"],
            df_m1_tg["end"],
            df_m1_tg["start"],
            df_m1_tg["strand"],
            df_m0_tg[f"{stat_to_average}_m0"],
            df_m1_tg[f"{stat_to_average}_m1"],
            df_m2_tg[f"{stat_to_average}_m2"],
            df_m3_tg[f"{stat_to_average}_m3"],
            df_m4_tg[f"{stat_to_average}_m4"],
            df_m5_tg[f"{stat_to_average}_m5"],
            df_m6_tg[f"{stat_to_average}_m6"],
            df_m7_tg[f"{stat_to_average}_m7"],
        ],
        axis=1,
    )

    # averaging over models
    stat_collected[f"{stat_to_average}"] = stat_collected[
        [
            f"{stat_to_average}_m{model_index}"
            for model_index in range(model_numbers)
        ]
    ].mean(axis=1)

    return stat_collected


def read_and_average_shuffling_exp(
    data_dir,
    not_shuffled_path="/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation/model_0.h5",
    stat_to_average="SCD",
    all_calculated_stats=["SCD", "INS-16", "INS-64"],
):
    # reading h5 files to dataframes
    df_no_shift = h5_to_df(
        not_shuffled_path, all_calculated_stats, average=False
    )

    df_n10k = h5_to_df(
        data_dir + "/model_0_shift_n10000.h5",
        all_calculated_stats,
        average=False,
    )
    df_n1k = h5_to_df(
        data_dir + "/model_0_shift_n1000.h5",
        all_calculated_stats,
        average=False,
    )
    df_n100 = h5_to_df(
        data_dir + "/model_0_shift_n100.h5",
        all_calculated_stats,
        average=False,
    )
    df_n10 = h5_to_df(
        data_dir + "/model_0_shift_n10.h5", all_calculated_stats, average=False
    )
    df_n1 = h5_to_df(
        data_dir + "/model_0_shift_n1.h5", all_calculated_stats, average=False
    )

    df_p10k = h5_to_df(
        data_dir + "/model_0_shift_p10000.h5",
        all_calculated_stats,
        average=False,
    )
    df_p1k = h5_to_df(
        data_dir + "/model_0_shift_p1000.h5",
        all_calculated_stats,
        average=False,
    )
    df_p100 = h5_to_df(
        data_dir + "/model_0_shift_p100.h5",
        all_calculated_stats,
        average=False,
    )
    df_p10 = h5_to_df(
        data_dir + "/model_0_shift_p10.h5", all_calculated_stats, average=False
    )
    df_p1 = h5_to_df(
        data_dir + "/model_0_shift_p1.h5", all_calculated_stats, average=False
    )

    # avergaing
    df_no_shift = average_stat_for_shift(
        df_no_shift, shift="no", model_index=0, head_index=1
    )

    df_n10k = average_stat_for_shift(
        df_n10k, shift="n10k", model_index=0, head_index=1
    )
    df_n1k = average_stat_for_shift(
        df_n1k, shift="n1k", model_index=0, head_index=1
    )
    df_n100 = average_stat_for_shift(
        df_n100, shift="n100", model_index=0, head_index=1
    )
    df_n10 = average_stat_for_shift(
        df_n10, shift="n10", model_index=0, head_index=1
    )
    df_n1 = average_stat_for_shift(
        df_n1, shift="n1", model_index=0, head_index=1
    )

    df_p10k = average_stat_for_shift(
        df_p10k, shift="p10k", model_index=0, head_index=1
    )
    df_p1k = average_stat_for_shift(
        df_p1k, shift="p1k", model_index=0, head_index=1
    )
    df_p100 = average_stat_for_shift(
        df_p100, shift="p100", model_index=0, head_index=1
    )
    df_p10 = average_stat_for_shift(
        df_p10, shift="p10", model_index=0, head_index=1
    )
    df_p1 = average_stat_for_shift(
        df_p1, shift="p1", model_index=0, head_index=1
    )

    stat_collected = pd.concat(
        [
            df_no_shift[f"{stat_to_average}_no"],
            df_n10k[f"{stat_to_average}_n10k"],
            df_n1k[f"{stat_to_average}_n1k"],
            df_n100[f"{stat_to_average}_n100"],
            df_n10[f"{stat_to_average}_n10"],
            df_n1[f"{stat_to_average}_n1"],
            df_p10k[f"{stat_to_average}_p10k"],
            df_p1k[f"{stat_to_average}_p1k"],
            df_p100[f"{stat_to_average}_p100"],
            df_p10[f"{stat_to_average}_p10"],
            df_p1[f"{stat_to_average}_p1"],
        ],
        axis=1,
    )

    return stat_collected


#########################









