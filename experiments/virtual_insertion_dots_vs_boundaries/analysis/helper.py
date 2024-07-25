import pandas as pd
from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import (average_stat_over_targets, average_stat_over_backgrounds)

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