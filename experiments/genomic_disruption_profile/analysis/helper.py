import pandas as pd
from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import average_stat_over_targets


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
