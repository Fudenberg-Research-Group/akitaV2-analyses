import numpy as np
import pandas as pd
import pysam
import scipy

from akita_utils.format_io import h5_to_df
from akita_utils.analysis_utils import (average_stat_over_targets,
                                        average_stat_over_backgrounds,
                                        average_stat_for_shift)


# HF5 reading functions

def read_and_average_virtual_exp(
    data_dir,
    stat_to_average="SCD",
    all_calculated_stats=["SCD", "INS-16", "INS-64"],
    model_numbers=8,
    head_index=1
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
    dfs = [h5_to_df(f"{data_dir}/model_{i}.h5", all_calculated_stats, average=False) for i in range(model_numbers)]
    print("averaging over targets")
    df_tg = [average_stat_over_targets(df, model_index=i, head_index=head_index, stat=stat_to_average) for i, df in enumerate(dfs)]

    print("averaging over backgrounds")
    df_tgbg = [average_stat_over_backgrounds(df, model_index=i, head_index=head_index, stat=stat_to_average) for i, df in enumerate(df_tg)]

    print(f"collecting data for {stat_to_average}")

    # Start with the first DataFrame and add the statistics columns
    stat_collected = df_tgbg[0][["chrom", "end", "start", "strand"]].copy()

    for i in range(model_numbers):
        stat_collected[f"{stat_to_average}_m{i}"] = df_tgbg[i][f"{stat_to_average}_m{i}"]

    # Ensure all columns are numeric and fill NaNs if needed
    stat_collected = stat_collected.fillna(0)  # You can use other strategies like forward/backward fill if needed

    # Averaging over models
    stat_collected[stat_to_average] = stat_collected[
        [f"{stat_to_average}_m{i}" for i in range(model_numbers)]
    ].mean(axis=1)

    return stat_collected


def read_and_average_genomic_exp(
    data_dir,
    stat_to_average="SCD",
    all_calculated_stats=["SCD", "INS-16", "INS-64"],
    model_numbers=8
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
    dfs = [h5_to_df(f"{data_dir}/model_{i}.h5", all_calculated_stats, average=False) for i in range(model_numbers)]

    print("averaging over targets")
    df_tg = [average_stat_over_targets(df, model_index=i, head_index=1, stat=stat_to_average) for i, df in enumerate(dfs)]

    print(f"collecting data for {stat_to_average}")
    # Collect data from all models
    data_frames = [df_tg[i][["chrom", "end", "start", "strand"]].assign(**{f"{stat_to_average}_m{i}": df_tg[i][f"{stat_to_average}_m{i}"]}) for i in range(model_numbers)]
    
    # Concatenate all DataFrames
    stat_collected = pd.concat(data_frames, axis=1)
    
    # Averaging over models
    stat_collected[stat_to_average] = stat_collected[
        [f"{stat_to_average}_m{i}" for i in range(model_numbers)]
    ].mean(axis=1)

    return stat_collected


def read_and_average_shifts(
    data_dir,
    not_shuffled_path="/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation/model_0.h5",
    stat_to_average="SCD",
    model_index=0,
    all_calculated_stats=["SCD", "INS-16", "INS-64"]
):
    """
    Reads data from h5 files for shuffled and non-shuffled genomic experiments, calculates the average statistic for each shift, and returns the concatenated result.

    Parameters:
    - data_dir (str): The directory path where the h5 files are stored.
    - not_shuffled_path (str): The path to the non-shuffled data file.
    - stat_to_average (str): The statistic to be averaged. Default is "SCD".
    - model_index (int): The index of the model. Default is 0.
    - all_calculated_stats (list of str): A list of all statistics calculated in the data. Default is ["SCD", "INS-16", "INS-64"].

    Returns:
    - DataFrame: A pandas DataFrame containing the concatenated and averaged statistics for each shift.
    """

    # Define shift types and shift values
    shifts = ['no', 'n1', 'n10', 'n100', 'n1000', 'n10000', 'p1', 'p10', 'p100', 'p1000', 'p10000']
    shift_paths = {
        'no': not_shuffled_path
    }
    shift_paths.update({
        shift: f"{data_dir}/model_{model_index}_shift_{shift}.h5" for shift in shifts[1:]
    })

    # Read and average statistics
    data_frames = {}
    for shift in shifts:
        df = h5_to_df(shift_paths[shift], all_calculated_stats, average=False)
        data_frames[shift] = average_stat_for_shift(df, shift=shift, model_index=model_index, head_index=1)

    # Collecting data
    stat_collected = pd.concat(
        [data_frames[shift][f"{stat_to_average}_{shift}"] for shift in shifts],
        axis=1
    )

    return stat_collected


def read_multi_model_shifts_data(data="/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/shifted_permutations/", num_models=4):
    """
    Reads and averages shifted data across multiple models and concatenates the results.
    
    Parameters:
    - data (str): Base directory for the model data.

    Returns:
    - DataFrame: A concatenated DataFrame with averaged statistics for each model.
    """

    # Initialize list to hold DataFrames for each model
    dfs = []

    # Process each model
    for model_index in range(num_models):  # Adjust the range if you have more models
        data_dir = f"{data}model_{model_index}"
        not_shuffled_path = f"/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation/model_{model_index}.h5"
        
        df = read_and_average_shifts(
            data_dir,
            not_shuffled_path=not_shuffled_path,
            stat_to_average="SCD",
            model_index=model_index
        )
        
        dfs.append(df)

    # Concatenate all DataFrames horizontally
    df_concat = pd.concat(dfs, axis=1)

    # Average columns across models
    cols = [col for col in df_concat.columns]
    for col in cols:
        df_concat[f"{col}_ave"] = df_concat.filter(like=col).mean(axis=1)

    return df_concat


def read_and_average_stat(data_dir, model_index, head_index, ignore_keys, stat, columns_to_keep, all_scores):
    """
    Reads and averages the given statistic from an H5 file for a specific model and head index.

    Parameters:
    - data_dir (str): The directory path where the H5 files are stored.
    - model_index (int): The index of the model.
    - head_index (int): The head index to use for averaging.
    - ignore_keys (list): List of keys to ignore while reading the H5 file.
    - stat (str): The statistic to be averaged.
    - columns_to_keep (list): List of columns to keep in the final DataFrame.

    Returns:
    - DataFrame: A pandas DataFrame containing the averaged statistic.
    """
    df = h5_to_df(data_dir, all_scores, ignore_keys=ignore_keys, average=False)
    df_tg = average_stat_over_targets(df, model_index=model_index, head_index=head_index, stat=stat)
    df_tgbg = average_stat_over_backgrounds(df_tg, model_index=model_index, head_index=head_index, stat=stat, columns_to_keep=columns_to_keep)
    df_tgbg = df_tgbg.rename(columns={f"{stat}_m{model_index}": f"{stat}_m{model_index}_D"})
    return df_tgbg


def read_average_dot_boundary(data_dir, model_index, head_index=1, ignore_keys=[], columns_to_keep=["chrom", "end", "start", "strand"], all_dots_scores=["SCD", "cross-score", "dot-score", "x-score"], all_boundary_scores=["SCD", "INS-16", "INS-64"]):
    """
    Reads and averages data for boundary and dot files, then combines them into a summary DataFrame.

    Parameters:
    - data_dir (str): The directory path where the H5 files are stored.
    - model_index (int): The index of the model.
    - head_index (int): The head index to use for averaging.
    - ignore_keys (list): List of keys to ignore while reading the H5 file.
    - columns_to_keep (list): List of columns to keep in the final DataFrame.

    Returns:
    - DataFrame: A pandas DataFrame containing the combined and averaged statistics.
    """
    # Process boundary data
    print("- Processing boundary scores")
    boundary_file = f"{data_dir}/model_{model_index}_boundary.h5"
    df_B = read_and_average_stat(boundary_file, model_index, head_index, ignore_keys, "SCD", columns_to_keep, all_boundary_scores)
    df_B = df_B.rename(columns={f"SCD_m{model_index}_D": f"SCD_m{model_index}_B"})

    # Process dot data for multiple stats
    print("- Processing dot scores")
    stats = ["SCD", "dot-score", "cross-score", "x-score"]
    dot_data = {}
    for stat in stats:
        dot_file = f"{data_dir}/model_{model_index}_dot.h5"
        df_dot = read_and_average_stat(dot_file, model_index, head_index, ignore_keys, stat, columns_to_keep, all_dots_scores)
        df_dot = df_dot.rename(columns={f"SCD_m{model_index}": f"SCD_m{model_index}_D"})
        dot_data[stat] = df_dot

    # Create summary DataFrame
    summary_df = df_B.drop(columns=[f"SCD_bg{x}" for x in range(10)])
    for stat in stats:
        summary_df[f"{stat}_m{model_index}"] = dot_data[stat][f"{stat}_m{model_index}_D"]
    summary_df = summary_df.rename(columns={f"SCD_m{model_index}": f"SCD_m{model_index}_D"})
    
    return summary_df


def summarize_average_models_dot_boundary(data_dir, models_number, head_index=1, ignore_keys=[], columns_to_keep=["chrom", "end", "start", "strand"], all_dots_scores=["SCD", "cross-score", "dot-score", "x-score"], all_boundary_scores=["SCD", "INS-16", "INS-64"]):

    print("Working on MODEL 0")
    df_summary = read_average_dot_boundary(data_dir, model_index=0, head_index=head_index,
                                          ignore_keys=ignore_keys, 
                                          columns_to_keep=columns_to_keep,
                                          all_dots_scores=all_dots_scores,
                                          all_boundary_scores=all_boundary_scores)
    
    for model_index in range(1, models_number):
        print(f"Working on MODEL {model_index}")
        model_df = read_average_dot_boundary(data_dir, model_index=model_index,
                                            head_index=head_index,
                                          ignore_keys=ignore_keys, 
                                          columns_to_keep=columns_to_keep,
                                          all_dots_scores=all_dots_scores,
                                          all_boundary_scores=all_boundary_scores)
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


# move dots here 


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
        df = h5_to_df(f"{data_dir}model_{i}.h5", ["SCD", "INS-16", "INS-64"], average=True)
        df = df.rename(columns={"SCD": f"SCD_m{i}", "INS-16": f"INS-16_m{i}", "INS-64": f"INS-64_m{i}"})
        if i == 0:
            df = df.drop(columns=["isBound", "name", "TF", "chipseq.score", "width"])
        dfs.append(df)

    # Merge all model DataFrames on the rownames column, with appropriate suffixes
    df_merged = dfs[0]
    for i in range(1, len(dfs)):
        df_merged = df_merged.merge(dfs[i], on="rownames", suffixes=('', f'_m{i}'))

    # Rename the rownames column to TFBS_cluster
    df_merged = df_merged.rename(columns={"rownames": "TFBS_cluster"})

    # Calculate the average for each statistic across models
    for stat in ["SCD", "INS-16", "INS-64"]:
        df_merged[stat] = df_merged[[f"{stat}_m{i}" for i in model_indices]].mean(axis=1)

    # Drop duplicates based on the TFBS_cluster column
    df_merged.drop_duplicates('TFBS_cluster', inplace=True)

    return df_merged


#############







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


def read_genomic_disruption_profile_data(genome_fragment, data_dir="/project/fudenber_735/akitaX1_analyses_data/genomic_disruption_profile/",
                                        mouse_target=0,
                                         human_target=1,
                                        mouse_head=1,
                                        human_head=0,
                                        columns_to_keep = ["chr", "start", "end", "perm_start", "perm_end"]):
    """
    Reads genomic disruption profile data for a given genome fragment, combining mouse and human data
    into a single dataframe.

    Parameters:
    -----------
    genome_fragment : str
        The name of the genome fragment to read data for.
    data_dir : str, optional
        The directory where the genomic disruption profile data is stored. Default is
        "/project/fudenber_735/akitaX1_analyses_data/genomic_disruption_profile/".
    mouse_target : int, optional
        The target index for the mouse data. Default is 0.
    human_target : int, optional
        The target index for the human data. Default is 1.
    mouse_head : int, optional
        The head index for the mouse data. Default is 1.
    human_head : int, optional
        The head index for the human data. Default is 0.
    columns_to_keep : list of str, optional
        The list of column names to keep in the resulting dataframe. Default is
        ["chr", "start", "end", "perm_start", "perm_end"].

    Returns:
    --------
    pd.DataFrame
        A dataframe containing the combined genomic disruption profile data for the specified genome
        fragment, with additional columns for mouse and human data and a "genomic_window_id" column
        representing unique genomic windows.
    """
    # preparing paths
    human_m0_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_human_m0.h5"
    human_m1_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_human_m1.h5"
    human_m2_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_human_m2.h5"
    human_m3_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_human_m3.h5"
    
    mouse_m0_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_mouse_m0.h5"
    mouse_m1_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_mouse_m1.h5"
    mouse_m2_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_mouse_m2.h5"
    mouse_m3_filename = f"{data_dir}/{genome_fragment}/{genome_fragment}_mouse_m3.h5"

    # preparing collective df
    mouse_m0_df = h5_to_df(mouse_m0_filename, stats=["SCD"], average=False, verbose=False)
    columns = columns_to_keep + [f"SCD_h{mouse_head}_m0_t{mouse_target}"]
    df = mouse_m0_df[columns]
    df = df.copy()
    df["genomic_window_id"] = pd.factorize(df['start'])[0]

    # reading the rest of mouse data
    mouse_m1_df = h5_to_df(mouse_m1_filename, stats=["SCD"], average=False, verbose=False)
    mouse_m2_df = h5_to_df(mouse_m2_filename, stats=["SCD"], average=False, verbose=False)
    mouse_m3_df = h5_to_df(mouse_m3_filename, stats=["SCD"], average=False, verbose=False)
    
    df[f"SCD_h{mouse_head}_m1_t{mouse_target}"] = mouse_m1_df[f"SCD_h{mouse_head}_m1_t{mouse_target}"]
    df[f"SCD_h{mouse_head}_m2_t{mouse_target}"] = mouse_m2_df[f"SCD_h{mouse_head}_m2_t{mouse_target}"]
    df[f"SCD_h{mouse_head}_m3_t{mouse_target}"] = mouse_m3_df[f"SCD_h{mouse_head}_m3_t{mouse_target}"]

    # reading human data and apending
    human_m0_df = h5_to_df(human_m0_filename, stats=["SCD"], average=False, verbose=False)
    human_m1_df = h5_to_df(human_m1_filename, stats=["SCD"], average=False, verbose=False)
    human_m2_df = h5_to_df(human_m2_filename, stats=["SCD"], average=False, verbose=False)
    human_m3_df = h5_to_df(human_m3_filename, stats=["SCD"], average=False, verbose=False)
    
    df[f"SCD_h{human_head}_m0_t{human_target}"] = human_m0_df[f"SCD_h{human_head}_m0_t{human_target}"]
    df[f"SCD_h{human_head}_m1_t{human_target}"] = human_m1_df[f"SCD_h{human_head}_m1_t{human_target}"]
    df[f"SCD_h{human_head}_m2_t{human_target}"] = human_m2_df[f"SCD_h{human_head}_m2_t{human_target}"]
    df[f"SCD_h{human_head}_m3_t{human_target}"] = human_m3_df[f"SCD_h{human_head}_m3_t{human_target}"]

    # computing averages for mouse and human columns
    df[f"mouse_avg_t{mouse_target}"] = df[[f"SCD_h{mouse_head}_m{i}_t{mouse_target}" for i in range(4)]].mean(axis=1)
    df[f"human_avg_t{human_target}"] = df[[f"SCD_h{human_head}_m{i}_t{human_target}" for i in range(4)]].mean(axis=1)
    
    return df

#########################









