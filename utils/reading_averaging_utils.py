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
    head_index=1
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
        df_m0, model_index=0, head_index=head_index, stat=stat_to_average
    )
    df_m1_tg = average_stat_over_targets(
        df_m1, model_index=1, head_index=head_index, stat=stat_to_average
    )
    df_m2_tg = average_stat_over_targets(
        df_m2, model_index=2, head_index=head_index, stat=stat_to_average
    )
    df_m3_tg = average_stat_over_targets(
        df_m3, model_index=3, head_index=head_index, stat=stat_to_average
    )
    df_m4_tg = average_stat_over_targets(
        df_m4, model_index=4, head_index=head_index, stat=stat_to_average
    )
    df_m5_tg = average_stat_over_targets(
        df_m5, model_index=5, head_index=head_index, stat=stat_to_average
    )
    df_m6_tg = average_stat_over_targets(
        df_m6, model_index=6, head_index=head_index, stat=stat_to_average
    )
    df_m7_tg = average_stat_over_targets(
        df_m7, model_index=7, head_index=head_index, stat=stat_to_average
    )

    print("averaging over backgrounds")
    df_m0_tgbg = average_stat_over_backgrounds(
        df_m0_tg, model_index=0, head_index=head_index, stat=stat_to_average
    )
    df_m1_tgbg = average_stat_over_backgrounds(
        df_m1_tg, model_index=1, head_index=head_index, stat=stat_to_average
    )
    df_m2_tgbg = average_stat_over_backgrounds(
        df_m2_tg, model_index=2, head_index=head_index, stat=stat_to_average
    )
    df_m3_tgbg = average_stat_over_backgrounds(
        df_m3_tg, model_index=3, head_index=head_index, stat=stat_to_average
    )
    df_m4_tgbg = average_stat_over_backgrounds(
        df_m4_tg, model_index=4, head_index=head_index, stat=stat_to_average
    )
    df_m5_tgbg = average_stat_over_backgrounds(
        df_m5_tg, model_index=5, head_index=head_index, stat=stat_to_average
    )
    df_m6_tgbg = average_stat_over_backgrounds(
        df_m6_tg, model_index=6, head_index=head_index, stat=stat_to_average
    )
    df_m7_tgbg = average_stat_over_backgrounds(
        df_m7_tg, model_index=7, head_index=head_index, stat=stat_to_average
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
    model_index=0,
    all_calculated_stats=["SCD", "INS-16", "INS-64"],
):
    # reading h5 files to dataframes
    df_no_shift = h5_to_df(
        not_shuffled_path, all_calculated_stats, average=False
    )

    df_n10k = h5_to_df(
        data_dir + f"/model_{model_index}_shift_n10000.h5",
        all_calculated_stats,
        average=False,
    )
    df_n1k = h5_to_df(
        data_dir + f"/model_{model_index}_shift_n1000.h5",
        all_calculated_stats,
        average=False,
    )
    df_n100 = h5_to_df(
        data_dir + f"/model_{model_index}_shift_n100.h5",
        all_calculated_stats,
        average=False,
    )
    df_n10 = h5_to_df(
        data_dir + f"/model_{model_index}_shift_n10.h5", all_calculated_stats, average=False
    )
    df_n1 = h5_to_df(
        data_dir + f"/model_{model_index}_shift_n1.h5", all_calculated_stats, average=False
    )

    df_p10k = h5_to_df(
        data_dir + f"/model_{model_index}_shift_p10000.h5",
        all_calculated_stats,
        average=False,
    )
    df_p1k = h5_to_df(
        data_dir + f"/model_{model_index}_shift_p1000.h5",
        all_calculated_stats,
        average=False,
    )
    df_p100 = h5_to_df(
        data_dir + f"/model_{model_index}_shift_p100.h5",
        all_calculated_stats,
        average=False,
    )
    df_p10 = h5_to_df(
        data_dir + f"/model_{model_index}_shift_p10.h5", all_calculated_stats, average=False
    )
    df_p1 = h5_to_df(
        data_dir + f"/model_{model_index}_shift_p1.h5", all_calculated_stats, average=False
    )

    # avergaing
    df_no_shift = average_stat_for_shift(
        df_no_shift, shift="no", model_index=model_index, head_index=1
    )

    df_n10k = average_stat_for_shift(
        df_n10k, shift="n10k", model_index=model_index, head_index=1
    )
    df_n1k = average_stat_for_shift(
        df_n1k, shift="n1k", model_index=model_index, head_index=1
    )
    df_n100 = average_stat_for_shift(
        df_n100, shift="n100", model_index=model_index, head_index=1
    )
    df_n10 = average_stat_for_shift(
        df_n10, shift="n10", model_index=model_index, head_index=1
    )
    df_n1 = average_stat_for_shift(
        df_n1, shift="n1", model_index=model_index, head_index=1
    )

    df_p10k = average_stat_for_shift(
        df_p10k, shift="p10k", model_index=model_index, head_index=1
    )
    df_p1k = average_stat_for_shift(
        df_p1k, shift="p1k", model_index=model_index, head_index=1
    )
    df_p100 = average_stat_for_shift(
        df_p100, shift="p100", model_index=model_index, head_index=1
    )
    df_p10 = average_stat_for_shift(
        df_p10, shift="p10", model_index=model_index, head_index=1
    )
    df_p1 = average_stat_for_shift(
        df_p1, shift="p1", model_index=model_index, head_index=1
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


def read_average_dot_boundary(data_dir, model_index, head_index=1, ignore_keys=[], columns_to_keep=["chrom", "end", "start", "strand"]):

    # reading and averaging data for model 0 only, boundaries
    df_B = h5_to_df(data_dir+f"/model_{model_index}_boundary.h5", ["SCD", "INS-16", "INS-64"], ignore_keys=ignore_keys, average=False) 
    df_B_tg = average_stat_over_targets(df_B, model_index=model_index, head_index=head_index, stat="SCD")
    df_B_tgbg = average_stat_over_backgrounds(df_B_tg, model_index=model_index, head_index=head_index, stat="SCD",
                                            columns_to_keep=columns_to_keep)

    df_B_tgbg = df_B_tgbg.rename(columns={f"SCD_m{model_index}": f"SCD_m{model_index}_B"})

    # reading and averaging data for model 0 only, dots
    df_D = h5_to_df(data_dir+f"/model_{model_index}_dot.h5", ["SCD", "dot-score", "cross-score", "x-score"], 
                       ignore_keys=ignore_keys, 
                       average=False)
    # SCD
    df_D_tg = average_stat_over_targets(df_D, model_index=model_index, head_index=head_index, stat="SCD")
    df_D_tgbg = average_stat_over_backgrounds(df_D_tg, model_index=model_index, head_index=head_index, stat="SCD",
                                            columns_to_keep=columns_to_keep)
    df_D_tgbg = df_D_tgbg.rename(columns={f"SCD_m{model_index}": f"SCD_m{model_index}_D"})

    # dot-score
    df_D_tg_dot = average_stat_over_targets(df_D, model_index=model_index, head_index=head_index, stat="dot-score")
    df_D_tgbg_dot = average_stat_over_backgrounds(df_D_tg_dot, model_index=model_index, head_index=head_index, stat="dot-score",
                                            columns_to_keep=columns_to_keep)

    # cross-score
    df_D_tg_cross = average_stat_over_targets(df_D, model_index=model_index, head_index=head_index, stat="cross-score")
    df_D_tgbg_cross = average_stat_over_backgrounds(df_D_tg_cross, model_index=model_index, head_index=head_index, stat="cross-score",
                                            columns_to_keep=columns_to_keep)

    # x-score
    df_D_tg_x = average_stat_over_targets(df_D, model_index=model_index, head_index=head_index, stat="x-score")
    df_D_tgbg_x = average_stat_over_backgrounds(df_D_tg_x, model_index=model_index, head_index=head_index, stat="x-score",
                                            columns_to_keep=columns_to_keep)

    # creating summary df
    summary_df = df_B_tgbg.drop(columns= [f"SCD_bg{x}" for x in range(10)])
    summary_df[f"SCD_m{model_index}_D"] = df_D_tgbg[f"SCD_m{model_index}_D"]
    summary_df[f"dot-score_m{model_index}"] = df_D_tgbg_dot[f"dot-score_m{model_index}"]
    summary_df[f"cross-score_m{model_index}"] = df_D_tgbg_cross[f"cross-score_m{model_index}"]
    summary_df[f"x-score_m{model_index}"] = df_D_tgbg_x[f"x-score_m{model_index}"]

    del df_D_tgbg
    del df_D_tgbg_dot
    del df_D_tgbg_cross
    del df_D_tgbg_x

    return summary_df


def summarize_average_models_dot_boundary(data_dir, models_number, head_index=1, ignore_keys=[], columns_to_keep=["chrom", "end", "start", "strand"]):
    df_summary = read_average_dot_boundary(data_dir, model_index=0, head_index=head_index,
                                          ignore_keys=ignore_keys, 
                                          columns_to_keep=columns_to_keep)
    
    for model_index in range(1, models_number):
        model_df = read_average_dot_boundary(data_dir, model_index=model_index,
                                            head_index=head_index,
                                          ignore_keys=ignore_keys, 
                                          columns_to_keep=columns_to_keep)
        df_summary = pd.concat([df_summary, 
                               model_df[f"SCD_m{model_index}_B"],
                               model_df[f"SCD_m{model_index}_D"],
                               model_df[f"dot-score_m{model_index}"],
                               model_df[f"cross-score_m{model_index}"],
                               model_df[f"x-score_m{model_index}"]], axis=1)
    
    df_summary["SCD_B"] = df_summary[[f"SCD_m{model_index}_B" for model_index in range(models_number)]].mean(axis=1)
    df_summary["SCD_D"] = df_summary[[f"SCD_m{model_index}_D" for model_index in range(models_number)]].mean(axis=1)
    df_summary["dot-score"] = df_summary[[f"dot-score_m{model_index}" for model_index in range(models_number)]].mean(axis=1)
    df_summary["cross-score"] = df_summary[[f"cross-score_m{model_index}" for model_index in range(models_number)]].mean(axis=1)
    df_summary["x-score"] = df_summary[[f"x-score_m{model_index}" for model_index in range(models_number)]].mean(axis=1)

    return df_summary


def read_disruption_smf_data(data_dir):
    
    # reading and averaging data for model 0 only
    df_m0 = h5_to_df(data_dir+"model_0.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m0 = df_m0.rename(columns={"SCD": "SCD_m0", "INS-16": "INS-16_m0", "INS-64": "INS-64_m0", "score": "PWM_score"})
    df_m0 = df_m0.drop(columns=["isBound", "name", "TF", "chipseq.score", "width"])

    # model 1
    df_m1 = h5_to_df(data_dir+"model_1.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m0["SCD_m1"] = df_m1["SCD"]
    df_m0["INS-16_m1"] = df_m1["INS-16"]
    df_m0["INS-64_m1"] = df_m1["INS-64"]

    # model 2
    df_m2 = h5_to_df(data_dir+"model_2.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m0["SCD_m2"] = df_m2["SCD"]
    df_m0["INS-16_m2"] = df_m2["INS-16"]
    df_m0["INS-64_m2"] = df_m2["INS-64"]

    # model 3
    df_m3 = h5_to_df(data_dir+"model_3.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m0["SCD_m3"] = df_m3["SCD"]
    df_m0["INS-16_m3"] = df_m3["INS-16"]
    df_m0["INS-64_m3"] = df_m3["INS-64"]
    
    df_m0.drop_duplicates('rownames',inplace=True)
    df_m0 = df_m0.rename(columns={"rownames": "TFBS_cluster"})
    
    # averaging
    df_m0["SCD"] = df_m0[["SCD_m0", "SCD_m1", "SCD_m2", "SCD_m3"]].mean(axis=1)
    df_m0["INS-16"] = df_m0[["INS-16_m0", "INS-16_m1", "INS-16_m2", "INS-16_m3"]].mean(axis=1)
    df_m0["INS-64"] = df_m0[["INS-64_m0", "INS-64_m1", "INS-64_m2", "INS-64_m3"]].mean(axis=1)

    return df_m0


def read_reverse_complement_disruption(disruption_dir, rc_disruption_dir):

    # model 0
    df_m0 = h5_to_df(disruption_dir+"model_0.h5", ["SCD", "INS-16", "INS-64"], average=True)
    df_rc_m0 = h5_to_df(rc_disruption_dir+"model_0_rc.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m0 = df_m0.rename(columns={"SCD": "SCD_m0", "INS-16": "INS-16_m0", "INS-64": "INS-64_m0"})
    df_rc_m0 = df_rc_m0.rename(columns={"SCD": "RC_SCD_m0", "INS-16": "RC_INS-16_m0", "INS-64": "RC_INS-64_m0"})
    df_m0[["RC_SCD_m0", "RC_INS-16_m0", "RC_INS-64_m0"]] = df_rc_m0[["RC_SCD_m0", "RC_INS-16_m0", "RC_INS-64_m0"]]

    # model 1
    df_m1 = h5_to_df(disruption_dir+"model_1.h5", ["SCD", "INS-16", "INS-64"], average=True)
    df_rc_m1 = h5_to_df(rc_disruption_dir+"model_1_rc.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m1 = df_m1.rename(columns={"SCD": "SCD_m1", "INS-16": "INS-16_m1", "INS-64": "INS-64_m1"})
    df_rc_m1 = df_rc_m1.rename(columns={"SCD": "RC_SCD_m1", "INS-16": "RC_INS-16_m1", "INS-64": "RC_INS-64_m1"})
    df_m0[["SCD_m1", "INS-16_m1", "INS-64_m1"]] = df_m1[["SCD_m1", "INS-16_m1", "INS-64_m1"]]
    df_m0[["RC_SCD_m1", "RC_INS-16_m1", "RC_INS-64_m1"]] = df_rc_m1[["RC_SCD_m1", "RC_INS-16_m1", "RC_INS-64_m1"]]

    # model 2
    df_m2 = h5_to_df(disruption_dir+"model_2.h5", ["SCD", "INS-16", "INS-64"], average=True)
    df_rc_m2 = h5_to_df(rc_disruption_dir+"model_2_rc.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m2 = df_m2.rename(columns={"SCD": "SCD_m2", "INS-16": "INS-16_m2", "INS-64": "INS-64_m2"})
    df_rc_m2 = df_rc_m2.rename(columns={"SCD": "RC_SCD_m2", "INS-16": "RC_INS-16_m2", "INS-64": "RC_INS-64_m2"})
    df_m0[["SCD_m2", "INS-16_m2", "INS-64_m2"]] = df_m2[["SCD_m2", "INS-16_m2", "INS-64_m2"]]
    df_m0[["RC_SCD_m2", "RC_INS-16_m2", "RC_INS-64_m2"]] = df_rc_m2[["RC_SCD_m2", "RC_INS-16_m2", "RC_INS-64_m2"]]

    # model 3
    df_m3 = h5_to_df(disruption_dir+"model_3.h5", ["SCD", "INS-16", "INS-64"], average=True)
    df_rc_m3 = h5_to_df(rc_disruption_dir+"model_3_rc.h5", ["SCD", "INS-16", "INS-64"], average=True) 
    df_m3 = df_m3.rename(columns={"SCD": "SCD_m3", "INS-16": "INS-16_m3", "INS-64": "INS-64_m3"})
    df_rc_m3 = df_rc_m3.rename(columns={"SCD": "RC_SCD_m3", "INS-16": "RC_INS-16_m3", "INS-64": "RC_INS-64_m3"})
    df_m0[["SCD_m3", "INS-16_m3", "INS-64_m3"]] = df_m3[["SCD_m3", "INS-16_m3", "INS-64_m3"]]
    df_m0[["RC_SCD_m3", "RC_INS-16_m3", "RC_INS-64_m3"]] = df_rc_m3[["RC_SCD_m3", "RC_INS-16_m3", "RC_INS-64_m3"]]

    # averaging
    df_m0["SCD"] = df_m0[["SCD_m0", "SCD_m1", "SCD_m2", "SCD_m3"]].mean(axis=1)
    df_m0["RC_SCD"] = df_m0[["RC_SCD_m0", "RC_SCD_m1", "RC_SCD_m2", "RC_SCD_m3"]].mean(axis=1)
    df_m0["INS-16"] = df_m0[["INS-16_m0", "INS-16_m1", "INS-16_m2", "INS-16_m3"]].mean(axis=1)
    df_m0["RC_INS-16"] = df_m0[["RC_INS-16_m0", "RC_INS-16_m1", "RC_INS-16_m2", "RC_INS-16_m3"]].mean(axis=1)
    df_m0["INS-64"] = df_m0[["INS-64_m0", "INS-64_m1", "INS-64_m2", "INS-64_m3"]].mean(axis=1)
    df_m0["RC_INS-64"] = df_m0[["RC_INS-64_m0", "RC_INS-64_m1", "RC_INS-64_m2", "RC_INS-64_m3"]].mean(axis=1)

    return df_m0


def read_shuffling_data(data_dir, stat_names=["SCD"]):
    # model 0
    df_m0 = h5_to_df(data_dir+"model_0.h5", stat_names, average=True) 
    df_m0 = df_m0.rename(columns={"SCD": "SCD_m0"})

    # model 1
    df_m1 = h5_to_df(data_dir+"model_1.h5", stat_names, average=True) 
    df_m1 = df_m1.rename(columns={"SCD": "SCD_m1"})
    df_m0["SCD_m1"] = df_m1["SCD_m1"]

    # model 2
    df_m2 = h5_to_df(data_dir+"model_2.h5", stat_names, average=True) 
    df_m2 = df_m2.rename(columns={"SCD": "SCD_m2"})
    df_m0["SCD_m2"] = df_m2["SCD_m2"]

    # model 3
    df_m3 = h5_to_df(data_dir+"model_3.h5", stat_names, average=True)
    df_m3 = df_m3.rename(columns={"SCD": "SCD_m3"})
    df_m0["SCD_m3"] = df_m3["SCD_m3"]

    # averaging
    df_m0["SCD"] = df_m0[["SCD_m0", "SCD_m1", "SCD_m2", "SCD_m3"]].mean(axis=1)

    return df_m0


def read_multi_model_shuffling_data(data="/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/shifted_permutations/"):

    # model 0
    model_index = 0
    data_dir = data+f"model_{model_index}"
    df_m0 = read_and_average_shuffling_exp(data_dir,
                                not_shuffled_path=f'/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation/model_{model_index}.h5',
                                    stat_to_average="SCD", model_index=model_index)  
    
    # model 1
    model_index = 1
    data_dir = data+f"model_{model_index}"
    df_m1 = read_and_average_shuffling_exp(data_dir,
                                not_shuffled_path=f'/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation/model_{model_index}.h5',
                                    stat_to_average="SCD", model_index=model_index)  

    # model 2
    model_index = 2
    data_dir = data+f"model_{model_index}"
    df_m2 = read_and_average_shuffling_exp(data_dir,
                                not_shuffled_path=f'/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation/model_{model_index}.h5',
                                    stat_to_average="SCD", model_index=model_index)  

    # model 3
    model_index = 3
    data_dir = data+f"model_{model_index}"
    df_m3 = read_and_average_shuffling_exp(data_dir,
                                not_shuffled_path=f'/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation/model_{model_index}.h5',
                                    stat_to_average="SCD", model_index=model_index)  

    df_concat = pd.concat([df_m0, df_m1, df_m2, df_m3], axis=1)  
    for col in df_concat.columns.unique():
        df_concat[f"{col}_ave"] = df_concat[[f"{col}"]].mean(axis=1)

    return df_concat


def read_multi_model_single_flanks_data(data_dir):
    
    # reading and averaging data for model 0, background 0 only
    df_m0_bg0 = h5_to_df(data_dir + "/model_0/" + "model_0_bg_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"]) 

    # preparing summary df
    df = average_stat_over_targets(df_m0_bg0, model_index=0, head_index=1, stat="SCD")
    df = df.rename(columns={"SCD_m0": "SCD_m0_bg0"})
    df["INS-16_m0_bg0"] = average_stat_over_targets(df_m0_bg0, model_index=0, head_index=1, stat="INS-16")["INS-16_m0"]
    df["INS-64_m0_bg0"] = average_stat_over_targets(df_m0_bg0, model_index=0, head_index=1, stat="INS-64")["INS-64_m0"]
    df = df.drop(columns=[f"SCD_h1_m0_t{x}" for x in range(6)] + [f"INS-16_h1_m0_t{x}" for x in range(6)] + [f"INS-64_h1_m0_t{x}" for x in range(6)] + ["background_index"])

    for model_index in range(4):
        print("Model: ", model_index)
        if model_index == 0:
            start_from = 1
        else:
            start_from = 0
    
        for bg_index in range(start_from, 10):
            print("\t - bg:", bg_index)
            df_bg = h5_to_df(data_dir + f"/model_{model_index}/" + f"model_{model_index}_bg_{bg_index}.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"]) 
            df_bg_tg = average_stat_over_targets(df_bg, model_index=model_index, head_index=1, stat="SCD")
            df_bg_tg = df_bg_tg.rename(columns={f"SCD_m{model_index}": f"SCD_m{model_index}_bg{bg_index}"})
            df[f"SCD_m{model_index}_bg{bg_index}"] = df_bg_tg[f"SCD_m{model_index}_bg{bg_index}"]
            df[f"INS-16_m{model_index}_bg{bg_index}"] = average_stat_over_targets(df_bg, model_index=model_index, head_index=1, stat="INS-16")[f"INS-16_m{model_index}"]
            df[f"INS-64_m{model_index}_bg{bg_index}"] = average_stat_over_targets(df_bg, model_index=model_index, head_index=1, stat="INS-64")[f"INS-64_m{model_index}"]

    # averaging
    df["SCD"] = df[[f"SCD_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)
    df["INS-16"] = df[[f"INS-16_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)
    df["INS-64"] = df[[f"INS-64_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)

    return df


def read_multi_model_double_flanks_data(data_dir):
    # model 0
    df_m0_0 = h5_to_df(data_dir + "model_0_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m0_1 = h5_to_df(data_dir + "model_0_1.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m0 = pd.concat([df_m0_0, df_m0_1])
    df_m0 = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat="SCD")
    df_m0["INS-16_m0"] = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat="INS-16")["INS-16_m0"]
    df_m0["INS-64_m0"] = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat="INS-64")["INS-64_m0"]
    df_m0 = df_m0.drop(columns=[f"SCD_h1_m0_t{target_ind}" for target_ind in range(6)] + [f"INS-16_h1_m0_t{target_ind}" for target_ind in range(6)] + [f"INS-64_h1_m0_t{target_ind}" for target_ind in range(6)])

    # model 1
    df_m1_0 = h5_to_df(data_dir + "model_1_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m1_1 = h5_to_df(data_dir + "model_1_1.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m1 = pd.concat([df_m1_0, df_m1_1])
    df_m0["SCD_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat="SCD")["SCD_m1"]
    df_m0["INS-16_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat="INS-16")["INS-16_m1"]
    df_m0["INS-64_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat="INS-64")["INS-64_m1"]
    del df_m1

    # model 2
    df_m2_0 = h5_to_df(data_dir + "model_2_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m2_1 = h5_to_df(data_dir + "model_2_1.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m2 = pd.concat([df_m2_0, df_m2_1])
    df_m0["SCD_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat="SCD")["SCD_m2"]
    df_m0["INS-16_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat="INS-16")["INS-16_m2"]
    df_m0["INS-64_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat="INS-64")["INS-64_m2"]
    del df_m2

    # model 3
    df_m3_0 = h5_to_df(data_dir + "model_3_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m3_1 = h5_to_df(data_dir + "model_3_1.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=["insertion_SCD", "disruption_SCD"])
    df_m3 = pd.concat([df_m3_0, df_m3_1])
    df_m0["SCD_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat="SCD")["SCD_m3"]
    df_m0["INS-16_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat="INS-16")["INS-16_m3"]
    df_m0["INS-64_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat="INS-64")["INS-64_m3"]
    del df_m3

    df = df_m0.groupby(["chrom", "start", "end", "strand", "orientation", "flank_bp"]).agg({"SCD_m0": "mean",
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
                                                                                  "INS-64_m3": "mean"}).reset_index()

    # averaging
    df["SCD"] = df[[f"SCD_m{model_ind}" for model_ind in range(4)]].mean(axis=1)
    df["INS-16"] = df[[f"INS-16_m{model_ind}" for model_ind in range(4)]].mean(axis=1)
    df["INS-64"] = df[[f"INS-64_m{model_ind}" for model_ind in range(4)]].mean(axis=1)
    return df


def read_multi_model_data(data_dir, keys_to_ignore = ["disruption_SCD_core",
                                                                     "disruption_SCD_flank",
                                                                     "insSCD_group_core",
                                                                     "insSCD_group_flank",
                                                                     "insertion_SCD_core",
                                                                     "insertion_SCD_flank"]):
    
    # reading and averaging data for model 0, background 0 only
    df_m0_bg0 = h5_to_df(data_dir + "/model_0/" + "model_0_bg0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=keys_to_ignore) 

    # preparing summary df
    df = average_stat_over_targets(df_m0_bg0, model_index=0, head_index=1, stat="SCD")
    df = df.rename(columns={"SCD_m0": "SCD_m0_bg0"})
    df["INS-16_m0_bg0"] = average_stat_over_targets(df_m0_bg0, model_index=0, head_index=1, stat="INS-16")["INS-16_m0"]
    df["INS-64_m0_bg0"] = average_stat_over_targets(df_m0_bg0, model_index=0, head_index=1, stat="INS-64")["INS-64_m0"]
    df = df.drop(columns=[f"SCD_h1_m0_t{x}" for x in range(6)] + [f"INS-16_h1_m0_t{x}" for x in range(6)] + [f"INS-64_h1_m0_t{x}" for x in range(6)] + ["background_index"])

    for model_index in range(4):
        print("Model: ", model_index)
        if model_index == 0:
            start_from = 1
        else:
            start_from = 0
    
        for bg_index in range(start_from, 10):
            print("\t - bg:", bg_index)
            df_bg = h5_to_df(data_dir + f"/model_{model_index}/" + f"model_{model_index}_bg{bg_index}.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=keys_to_ignore) 
            df_bg_tg = average_stat_over_targets(df_bg, model_index=model_index, head_index=1, stat="SCD")
            df_bg_tg = df_bg_tg.rename(columns={f"SCD_m{model_index}": f"SCD_m{model_index}_bg{bg_index}"})
            df[f"SCD_m{model_index}_bg{bg_index}"] = df_bg_tg[f"SCD_m{model_index}_bg{bg_index}"]
            df[f"INS-16_m{model_index}_bg{bg_index}"] = average_stat_over_targets(df_bg, model_index=model_index, head_index=1, stat="INS-16")[f"INS-16_m{model_index}"]
            df[f"INS-64_m{model_index}_bg{bg_index}"] = average_stat_over_targets(df_bg, model_index=model_index, head_index=1, stat="INS-64")[f"INS-64_m{model_index}"]

    # averaging
    df["SCD"] = df[[f"SCD_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)
    df["INS-16"] = df[[f"INS-16_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)
    df["INS-64"] = df[[f"INS-64_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)

    return df


def read_multi_model_single_mutagenesis_data(data_dir, keys_to_ignore = ["disruption_SCD", "insertion_SCD"]):
    
    # model 0
    df_m0 = h5_to_df(data_dir+"/model_0.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=keys_to_ignore) 
    df = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat="SCD")
    df["INS-16_m0"] = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat="INS-16")["INS-16_m0"]
    df["INS-64_m0"] = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat="INS-64")["INS-64_m0"]
    df = df.drop(columns=[f"SCD_h1_m0_t{x}" for x in range(6)] + [f"INS-16_h1_m0_t{x}" for x in range(6)] + [f"INS-64_h1_m0_t{x}" for x in range(6)] + ["background_index"])

    # model 1
    df_m1 = h5_to_df(data_dir+"/model_1.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=keys_to_ignore) 
    df["SCD_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat="SCD")["SCD_m1"]
    df["INS-16_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat="INS-16")["INS-16_m1"]
    df["INS-64_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat="INS-64")["INS-64_m1"]

    # model 2
    df_m2 = h5_to_df(data_dir+"/model_2.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=keys_to_ignore) 
    df["SCD_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat="SCD")["SCD_m2"]
    df["INS-16_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat="INS-16")["INS-16_m2"]
    df["INS-64_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat="INS-64")["INS-64_m2"]

    # model 3
    df_m3 = h5_to_df(data_dir+"/model_3.h5", ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=keys_to_ignore) 
    df["SCD_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat="SCD")["SCD_m3"]
    df["INS-16_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat="INS-16")["INS-16_m3"]
    df["INS-64_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat="INS-64")["INS-64_m3"]

    # averaging
    df["SCD"] = df[[f"SCD_m{model_ind}" for model_ind in range(4)]].mean(axis=1)
    df["INS-16"] = df[[f"INS-16_m{model_ind}" for model_ind in range(4)]].mean(axis=1)
    df["INS-64"] = df[[f"INS-64_m{model_ind}" for model_ind in range(4)]].mean(axis=1)

    return df


def read_and_average_insulation_offset_data(data_dir, keys_to_ignore=["insertion_SCD", "disruption_SCD"]):
    """
    Reads insulation and offset data from a series of HDF5 files within a specified directory, 
    then calculates and appends averages across multiple models for selected statistics.
    
    This function operates by loading data from four model-specific HDF5 files named 'model_0.h5' 
    through 'model_3.h5'. It ignores specified keys in the data, then averages specified statistics
    (OFF-16, OFF-64, and OFF-128) over targets within each model. These averages are calculated
    separately for each model and appended to a DataFrame. Finally, the function computes the overall
    mean of each statistic across all models and returns a DataFrame with these values alongside 
    the model-specific averages.

    Parameters:
    - data_dir (str): The directory path where the model HDF5 files are located. It is expected that
                      the files are named sequentially as 'model_0.h5' to 'model_3.h5'.
    - keys_to_ignore (list of str, optional): A list of keys to ignore when reading data from the HDF5
                      files. Defaults to ["insertion_SCD", "disruption_SCD"].

    Returns:
    - pandas.DataFrame: A DataFrame containing the averaged statistics (OFF-16, OFF-64, and OFF-128) for
                        each model (as separate columns) and the overall averages of these statistics 
                        across all models (as separate columns).
    """
    
    df_m0 = h5_to_df(data_dir+"model_0.h5", ["SCD", "OFF-16", "OFF-64", "OFF-128"], ignore_keys=keys_to_ignore, average=False)
    df_m1 = h5_to_df(data_dir+"model_1.h5", ["SCD", "OFF-16", "OFF-64", "OFF-128"], ignore_keys=keys_to_ignore, average=False)
    df_m2 = h5_to_df(data_dir+"model_2.h5", ["SCD", "OFF-16", "OFF-64", "OFF-128"], ignore_keys=keys_to_ignore, average=False)
    df_m3 = h5_to_df(data_dir+"model_3.h5", ["SCD", "OFF-16", "OFF-64", "OFF-128"], ignore_keys=keys_to_ignore, average=False)

    df = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat='OFF-64')
    df["OFF-16_m0"] = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat='OFF-16')["OFF-16_m0"]
    df["OFF-128_m0"] = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat='OFF-128')["OFF-128_m0"]

    df["OFF-16_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat='OFF-16')["OFF-16_m1"]
    df["OFF-64_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat='OFF-64')["OFF-64_m1"]
    df["OFF-128_m1"] = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat='OFF-128')["OFF-128_m1"]

    df["OFF-16_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat='OFF-16')["OFF-16_m2"]
    df["OFF-64_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat='OFF-64')["OFF-64_m2"]
    df["OFF-128_m2"] = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat='OFF-128')["OFF-128_m2"]

    df["OFF-16_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat='OFF-16')["OFF-16_m3"]
    df["OFF-64_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat='OFF-64')["OFF-64_m3"]
    df["OFF-128_m3"] = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat='OFF-128')["OFF-128_m3"]

    df = df.copy()
    df["OFF-16"] = df[[f"OFF-16_m{i}" for i in range(4)]].mean(axis=1)
    df["OFF-64"] = df[[f"OFF-64_m{i}" for i in range(4)]].mean(axis=1)
    df["OFF-128"] = df[[f"OFF-128_m{i}" for i in range(4)]].mean(axis=1)

    return df


def summarize_dot_anchors_data(data_dir, head_index=1, columns_to_keep=["chrom", "end", "start", "strand"], ignore_keys=[]):

    # BOUNDARY DATA
    print("- processing boundary data from model 0")
    df_B_bg04 = h5_to_df(data_dir+f"/model_0_boundary_bg04.h5", ["SCD", "INS-16", "INS-64"], ignore_keys=ignore_keys, average=False) 
    df_B_bg59 = h5_to_df(data_dir+f"/model_0_boundary_bg59.h5", ["SCD", "INS-16", "INS-64"], ignore_keys=ignore_keys, average=False) 
    df_B = pd.concat([df_B_bg04, df_B_bg59]).reset_index(drop=True)

    # averaging over targets
    df_B_ave = average_stat_over_targets(df_B, model_index=0, head_index=1, stat="SCD")
    df_B_ave["INS-16_m0"] = average_stat_over_targets(df_B, model_index=0, head_index=head_index, stat="INS-16")["INS-16_m0"]
    df_B_ave["INS-64_m0"] = average_stat_over_targets(df_B, model_index=0, head_index=head_index, stat="INS-64")["INS-64_m0"]

    # averaging over backgrounds
    summary_df = average_stat_over_backgrounds(df_B_ave, model_index=0, head_index=head_index, stat="SCD",columns_to_keep=columns_to_keep)[["chrom", "end", "start", "strand", "SCD_m0"]]
    summary_df = summary_df.rename(columns={"SCD_m0": "SCD_B_m0"})
    summary_df["INS-16_B_m0"] = average_stat_over_backgrounds(df_B_ave, model_index=0, head_index=head_index, stat="INS-16",columns_to_keep=columns_to_keep)["INS-16_m0"]
    summary_df["INS-64_B_m0"] = average_stat_over_backgrounds(df_B_ave, model_index=0, head_index=head_index, stat="INS-64",columns_to_keep=columns_to_keep)["INS-64_m0"]

    # adding data from the models 1,2,3
    for model_index in range(1,4):
        print("- processing boundary data from model", model_index)
        data_bg04 = h5_to_df(data_dir+f"/model_{model_index}_boundary_bg04.h5", ["SCD", "INS-16", "INS-64"], ignore_keys=ignore_keys, average=False) 
        data_bg59 = h5_to_df(data_dir+f"/model_{model_index}_boundary_bg59.h5", ["SCD", "INS-16", "INS-64"], ignore_keys=ignore_keys, average=False) 
        data = pd.concat([data_bg04, data_bg59]).reset_index(drop=True)
        
        # averaging over targets
        df_B_ave[f"SCD_m{model_index}"] = average_stat_over_targets(data, model_index=model_index, head_index=head_index, stat="SCD")[f"SCD_m{model_index}"]
        df_B_ave[f"INS-16_m{model_index}"] = average_stat_over_targets(data, model_index=model_index, head_index=head_index, stat="INS-16")[f"INS-16_m{model_index}"]
        df_B_ave[f"INS-64_m{model_index}"] = average_stat_over_targets(data, model_index=model_index, head_index=head_index, stat="INS-64")[f"INS-64_m{model_index}"]
    
        # averaging over backgrounds
        summary_df[f"SCD_B_m{model_index}"] = average_stat_over_backgrounds(df_B_ave, model_index=model_index, head_index=head_index, stat="SCD",columns_to_keep=columns_to_keep)[f"SCD_m{model_index}"]
        summary_df[f"INS-16_B_m{model_index}"] = average_stat_over_backgrounds(df_B_ave, model_index=model_index, head_index=head_index, stat="INS-16",columns_to_keep=columns_to_keep)[f"INS-16_m{model_index}"]
        summary_df[f"INS-64_B_m{model_index}"] = average_stat_over_backgrounds(df_B_ave, model_index=model_index, head_index=head_index, stat="INS-64",columns_to_keep=columns_to_keep)[f"INS-64_m{model_index}"]


    # DOT DATA
    print("- processing dot data from model 0")
    df_D_bg04 = h5_to_df(data_dir+f"/model_0_dot_bg04.h5", ["SCD", "dot-score", "cross-score", "x-score"], 
                   ignore_keys=ignore_keys, 
                       average=False)
    df_D_bg59 = h5_to_df(data_dir+f"/model_0_dot_bg59.h5", ["SCD", "dot-score", "cross-score", "x-score"], 
                   ignore_keys=ignore_keys, 
                       average=False)
    df_D = pd.concat([df_D_bg04, df_D_bg59]).reset_index(drop=True)

    # averaging over targets
    df_D_ave = average_stat_over_targets(df_D, model_index=0, head_index=head_index, stat="SCD")
    df_D_ave["dot-score_m0"] = average_stat_over_targets(df_D, model_index=0, head_index=head_index, stat="dot-score")["dot-score_m0"]
    df_D_ave["cross-score_m0"] = average_stat_over_targets(df_D, model_index=0, head_index=head_index, stat="cross-score")["cross-score_m0"]
    df_D_ave["x-score_m0"] = average_stat_over_targets(df_D, model_index=0, head_index=head_index, stat="x-score")["x-score_m0"]

    # averaging over backgrounds
    summary_df["SCD_D_m0"] = average_stat_over_backgrounds(df_D_ave, model_index=0, head_index=head_index, stat="SCD", columns_to_keep=columns_to_keep)[["SCD_m0"]]
    summary_df["dot-score_D_m0"] = average_stat_over_backgrounds(df_D_ave, model_index=0, head_index=head_index, stat="dot-score", columns_to_keep=columns_to_keep)["dot-score_m0"]
    summary_df["cross-score_D_m0"] = average_stat_over_backgrounds(df_D_ave, model_index=0, head_index=head_index, stat="cross-score", columns_to_keep=columns_to_keep)["cross-score_m0"]
    summary_df["x-score_D_m0"] = average_stat_over_backgrounds(df_D_ave, model_index=0, head_index=head_index, stat="x-score", columns_to_keep=columns_to_keep)["x-score_m0"]

    # adding data from the models 1,2,3
    for model_index in range(1,4):
        print("- processing dot data from model", model_index)
        data_bg04 = h5_to_df(data_dir+f"/model_{model_index}_dot_bg04.h5", ["SCD", "dot-score", "cross-score", "x-score"], average=False) 
        data_bg59 = h5_to_df(data_dir+f"/model_{model_index}_dot_bg59.h5", ["SCD", "dot-score", "cross-score", "x-score"], average=False) 
        data = pd.concat([data_bg04, data_bg59]).reset_index(drop=True)
        
        # averaging over targets
        df_D_ave[f"SCD_m{model_index}"] = average_stat_over_targets(data, model_index=model_index, head_index=head_index, stat="SCD")[f"SCD_m{model_index}"]
        df_D_ave[f"dot-score_m{model_index}"] = average_stat_over_targets(data, model_index=model_index, head_index=head_index, stat="dot-score")[f"dot-score_m{model_index}"]
        df_D_ave[f"cross-score_m{model_index}"] = average_stat_over_targets(data, model_index=model_index, head_index=head_index, stat="cross-score")[f"cross-score_m{model_index}"]
        df_D_ave[f"x-score_m{model_index}"] = average_stat_over_targets(data, model_index=model_index, head_index=head_index, stat="x-score")[f"x-score_m{model_index}"]
    
        # averaging over backgrounds
        summary_df[[f"SCD_D_m{model_index}", "DETECTION_SCALE", "FDR"]] = average_stat_over_backgrounds(df_D_ave, model_index=model_index, head_index=head_index, stat="SCD",columns_to_keep=columns_to_keep+["DETECTION_SCALE", "FDR"])[[f"SCD_m{model_index}", "DETECTION_SCALE", "FDR"]]
        summary_df[f"dot-score_D_m{model_index}"] = average_stat_over_backgrounds(df_D_ave, model_index=model_index, head_index=head_index, stat="dot-score",columns_to_keep=columns_to_keep)[f"dot-score_m{model_index}"]
        summary_df[f"cross-score_D_m{model_index}"] = average_stat_over_backgrounds(df_D_ave, model_index=model_index, head_index=head_index, stat="cross-score",columns_to_keep=columns_to_keep)[f"cross-score_m{model_index}"]
        summary_df[f"x-score_D_m{model_index}"] = average_stat_over_backgrounds(df_D_ave, model_index=model_index, head_index=head_index, stat="x-score",columns_to_keep=columns_to_keep)[f"x-score_m{model_index}"]


    # averaging over models
    summary_df["SCD_B"] = summary_df[[f"SCD_B_m{i}" for i in range(4)]].mean(axis=1)
    summary_df["INS-16_B"] = summary_df[[f"INS-16_B_m{i}" for i in range(4)]].mean(axis=1)
    summary_df["INS-64_B"] = summary_df[[f"INS-64_B_m{i}" for i in range(4)]].mean(axis=1)

    summary_df["SCD_D"] = summary_df[[f"SCD_D_m{i}" for i in range(4)]].mean(axis=1)
    summary_df["dot-score_D"] = summary_df[[f"dot-score_D_m{i}" for i in range(4)]].mean(axis=1)
    summary_df["cross-score_D"] = summary_df[[f"cross-score_D_m{i}" for i in range(4)]].mean(axis=1)
    summary_df["x-score_D"] = summary_df[[f"x-score_D_m{i}" for i in range(4)]].mean(axis=1)

    return summary_df

#########################









