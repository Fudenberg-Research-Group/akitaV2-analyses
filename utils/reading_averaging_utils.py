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


def summarize_dot_anchors_data(data_dir, head_index=1, columns_to_keep=["chrom", "end", "start", "strand"], ignore_keys=[]):
    def process_data(model_index, data_type, stats, include_extra_columns=False):
        print(f"- processing {data_type} data from model {model_index}")
        df_bg04 = h5_to_df(f"{data_dir}/model_{model_index}_{data_type}_bg04.h5", stats, ignore_keys=ignore_keys, average=False)
        df_bg59 = h5_to_df(f"{data_dir}/model_{model_index}_{data_type}_bg59.h5", stats, ignore_keys=ignore_keys, average=False)
        df = pd.concat([df_bg04, df_bg59]).reset_index(drop=True)
        
        df_ave = average_stat_over_targets(df, model_index=model_index, head_index=head_index, stat=stats[0])
        for stat in stats[1:]:
            df_ave[f"{stat}_m{model_index}"] = average_stat_over_targets(df, model_index=model_index, head_index=head_index, stat=stat)[f"{stat}_m{model_index}"]

        if include_extra_columns:
            extra_columns = ["DETECTION_SCALE", "FDR"]
            for col in extra_columns:
                df_ave[col] = df[col]

        return df_ave

    def add_data_to_summary(df_ave, model_index, data_type, stats):
        for stat in stats:
            summary_df[f"{stat}_{data_type}_m{model_index}"] = average_stat_over_backgrounds(df_ave, model_index=model_index, head_index=head_index, stat=stat, columns_to_keep=columns_to_keep)[f"{stat}_m{model_index}"]

    # Boundary Data
    boundary_stats = ["SCD", "INS-16", "INS-64"]
    df_B_ave = process_data(0, "boundary", boundary_stats)
    summary_df = average_stat_over_backgrounds(df_B_ave, model_index=0, head_index=head_index, stat="SCD", columns_to_keep=columns_to_keep)[["chrom", "end", "start", "strand", "SCD_m0"]]
    summary_df = summary_df.rename(columns={"SCD_m0": "SCD_B_m0"})
    for stat in boundary_stats[1:]:
        summary_df[f"{stat}_B_m0"] = average_stat_over_backgrounds(df_B_ave, model_index=0, head_index=head_index, stat=stat, columns_to_keep=columns_to_keep)[f"{stat}_m0"]

    for model_index in range(1, 4):
        df_B_ave = process_data(model_index, "boundary", boundary_stats)
        add_data_to_summary(df_B_ave, model_index, "B", boundary_stats)

    # Dot Data
    dot_stats = ["SCD", "dot-score", "cross-score", "x-score"]
    df_D_ave = process_data(0, "dot", dot_stats, include_extra_columns=True)
    for stat in dot_stats:
        summary_df[f"{stat}_D_m0"] = average_stat_over_backgrounds(df_D_ave, model_index=0, head_index=head_index, stat=stat, columns_to_keep=columns_to_keep)[f"{stat}_m0"]

    summary_df["DETECTION_SCALE"] = df_D_ave["DETECTION_SCALE"]
    summary_df["FDR"] = df_D_ave["FDR"]

    for model_index in range(1, 4):
        df_D_ave = process_data(model_index, "dot", dot_stats)
        add_data_to_summary(df_D_ave, model_index, "D", dot_stats)

    # Averaging over models
    for stat in ["SCD_B", "INS-16_B", "INS-64_B", "SCD_D", "dot-score_D", "cross-score_D", "x-score_D"]:
        summary_df[stat] = summary_df[[f"{stat}_m{i}" for i in range(4)]].mean(axis=1)

    return summary_df


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


def read_genomic_disruption_profile_data(genome_fragment, data_dir="/project/fudenber_735/akitaX1_analyses_data/genomic_disruption_profile/",
                                         mouse_target=0, human_target=1, mouse_head=1, human_head=0,
                                         columns_to_keep=["chr", "start", "end", "perm_start", "perm_end"]):
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
    def load_and_merge_data(filenames, prefix, head, target):
        dfs = [h5_to_df(f, stats=["SCD"], average=False, verbose=False) for f in filenames]
        merged_df = dfs[0][columns_to_keep + [f"SCD_h{head}_m0_t{target}"]].copy()
        for i in range(1, 4):
            merged_df[f"SCD_h{head}_m{i}_t{target}"] = dfs[i][f"SCD_h{head}_m{i}_t{target}"]
        return merged_df

    mouse_files = [f"{data_dir}/{genome_fragment}/{genome_fragment}_mouse_m{i}.h5" for i in range(4)]
    human_files = [f"{data_dir}/{genome_fragment}/{genome_fragment}_human_m{i}.h5" for i in range(4)]
    
    df = load_and_merge_data(mouse_files, "mouse", mouse_head, mouse_target)
    df["genomic_window_id"] = pd.factorize(df['start'])[0]

    human_df = load_and_merge_data(human_files, "human", human_head, human_target)
    for i in range(4):
        df[f"SCD_h{human_head}_m{i}_t{human_target}"] = human_df[f"SCD_h{human_head}_m{i}_t{human_target}"]

    df[f"mouse_avg_t{mouse_target}"] = df[[f"SCD_h{mouse_head}_m{i}_t{mouse_target}" for i in range(4)]].mean(axis=1)
    df[f"human_avg_t{human_target}"] = df[[f"SCD_h{human_head}_m{i}_t{human_target}" for i in range(4)]].mean(axis=1)

    return df


# flank analysis

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


# flank-core compatibility & pairwise mutagenesis

def read_multi_model_compatibility_data(data_dir, keys_to_ignore=None):
    """
    Reads and processes data from multiple models and backgrounds to compute statistical summaries.

    This function processes HDF5 files from multiple models and background indices, calculates statistical
    summaries for different target categories, and returns a consolidated DataFrame with average statistics.

    Parameters:
    ----------
    data_dir : str
        Directory path where the HDF5 files for each model and background are located.
    
    keys_to_ignore : list of str, optional
        List of keys to ignore when reading HDF5 files. Default is a predefined list of keys related to disruptions
        and insertions.

    Returns:
    -------
    pandas.DataFrame
        A DataFrame containing the following columns:
        - `SCD_m0_bg0`, `INS-16_m0_bg0`, `INS-64_m0_bg0`: Summary statistics for model 0, background 0.
        - `SCD_m{model_index}_bg{bg_index}`, `INS-16_m{model_index}_bg{bg_index}`, `INS-64_m{model_index}_bg{bg_index}`:
          Summary statistics for each model and background combination.
        - `SCD`, `INS-16`, `INS-64`: Averaged statistics across all models and backgrounds.
    """
    if keys_to_ignore is None:
        keys_to_ignore = ["disruption_SCD_core", "disruption_SCD_flank", "insSCD_group_core",
                          "insSCD_group_flank", "insertion_SCD_core", "insertion_SCD_flank"]
    
    def process_data(model_index, bg_index):
        df_bg = h5_to_df(f"{data_dir}/model_{model_index}/model_{model_index}_bg{bg_index}.h5", 
                         ["SCD", "INS-16", "INS-64"], average=False, ignore_keys=keys_to_ignore)
        return {stat: average_stat_over_targets(df_bg, model_index=model_index, head_index=1, stat=stat)
                .rename(columns={f"{stat}_m{model_index}": f"{stat}_m{model_index}_bg{bg_index}"})
                for stat in ["SCD", "INS-16", "INS-64"]}

    # Initialize DataFrame with model 0, background 0 data
    df = h5_to_df(f"{data_dir}/model_0/model_0_bg0.h5", ["SCD", "INS-16", "INS-64"], 
                  average=False, ignore_keys=keys_to_ignore)
    df = average_stat_over_targets(df, model_index=0, head_index=1, stat="SCD")
    df = df.rename(columns={"SCD_m0": "SCD_m0_bg0"})
    df["INS-16_m0_bg0"] = average_stat_over_targets(df, model_index=0, head_index=1, stat="INS-16")["INS-16_m0"]
    df["INS-64_m0_bg0"] = average_stat_over_targets(df, model_index=0, head_index=1, stat="INS-64")["INS-64_m0"]
    df = df.drop(columns=[f"SCD_h1_m0_t{x}" for x in range(6)] +
                      [f"INS-16_h1_m0_t{x}" for x in range(6)] +
                      [f"INS-64_h1_m0_t{x}" for x in range(6)] + ["background_index"])

    # Process and add data for all models and backgrounds
    for model_index in range(4):
        start_from = 1 if model_index == 0 else 0
        for bg_index in range(start_from, 10):
            stats = process_data(model_index, bg_index)
            for stat, df_stat in stats.items():
                df[f"{stat}_m{model_index}_bg{bg_index}"] = df_stat[f"{stat}_m{model_index}_bg{bg_index}"]

    # Averaging
    for stat in ["SCD", "INS-16", "INS-64"]:
        df[stat] = df[[f"{stat}_m{model_ind}_bg{bg_ind}" for model_ind in range(4) for bg_ind in range(10)]].mean(axis=1)

    return df


# single mutagenesis

def read_multi_model_single_mutagenesis_data(data_dir, keys_to_ignore=["disruption_SCD", "insertion_SCD"]):
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
        model_df = h5_to_df(f"{data_dir}/model_{model_index}.h5", stats, average=False, ignore_keys=keys_to_ignore)
        if model_index == 0:
            df = model_df[["chrom", "start", "end", "strand", "disruption_SCD", "insertion_SCD", "mutated_nucleotide", "original_nucleotide", "position"]].copy()
        for stat in stats:
            df_stat = average_stat_over_targets(model_df, model_index=model_index, head_index=1, stat=stat)
            df[f"{stat}_m{model_index}"] = df_stat[f"{stat}_m{model_index}"]
    
    for stat in stats:
        df[stat] = df[[f"{stat}_m{model_index}" for model_index in range(4)]].mean(axis=1)
    
    return df


# insulation offset

def read_and_average_insulation_offset_data(data_dir, keys_to_ignore=["insertion_SCD", "disruption_SCD"]):
    """
    Reads insulation and offset data from a series of HDF5 files within a specified directory,
    then calculates and appends averages across multiple models for selected statistics.
    
    Parameters:
    - data_dir (str): Directory path where the model HDF5 files are located. Files are named 'model_0.h5' to 'model_3.h5'.
    - keys_to_ignore (list of str, optional): List of keys to ignore when reading data. Defaults to ["insertion_SCD", "disruption_SCD"].

    Returns:
    - pandas.DataFrame: DataFrame with averaged statistics (OFF-16, OFF-64, OFF-128) for each model and overall averages.
    """
    
    models = range(4)  # Model indices
    stats = ['OFF-16', 'OFF-64', 'OFF-128']  # Statistics to be averaged
    
    # Read data for each model into a dictionary of DataFrames
    dfs = {model: h5_to_df(f"{data_dir}model_{model}.h5", ["SCD"] + stats, ignore_keys=keys_to_ignore, average=False) for model in models}
    
    # Initialize DataFrame to store results
    # df = pd.DataFrame()
    df = dfs[0][["chrom", "start", "end", "strand", "orientation"]].copy()
    
    # Calculate averages for each statistic across all models
    for stat in stats:
        for model in models:
            df[f"{stat}_m{model}"] = average_stat_over_targets(dfs[model], model_index=model, head_index=1, stat=stat)[f"{stat}_m{model}"]
        
        # Compute overall average for each statistic
        df[stat] = df[[f"{stat}_m{model}" for model in models]].mean(axis=1)

    return df

