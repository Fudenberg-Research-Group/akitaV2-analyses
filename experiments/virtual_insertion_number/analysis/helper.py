from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import average_stat_over_targets


def read_multi_model_number_data(data_dir, keys_to_ignore=None):
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
        keys_to_ignore = [
            "disruption_SCD_core",
            "disruption_SCD_flank",
            "insSCD_group_core",
            "insSCD_group_flank",
            "insertion_SCD_core",
            "insertion_SCD_flank",
        ]

    def process_data(model_index, bg_index):
        df_bg = h5_to_df(
            f"{data_dir}/model_{model_index}/model_{model_index}_bg{bg_index}.h5",
            ["SCD", "INS-16", "INS-64"],
            average=False,
            ignore_keys=keys_to_ignore,
        )
        return {
            stat: average_stat_over_targets(
                df_bg, model_index=model_index, head_index=1, stat=stat
            ).rename(
                columns={
                    f"{stat}_m{model_index}": f"{stat}_m{model_index}_bg{bg_index}"
                }
            )
            for stat in ["SCD", "INS-16", "INS-64"]
        }

    # Initialize DataFrame with model 0, background 0 data
    df = h5_to_df(
        f"{data_dir}/model_0/model_0_bg0.h5",
        ["SCD", "INS-16", "INS-64"],
        average=False,
        ignore_keys=keys_to_ignore,
    )
    df = average_stat_over_targets(df, model_index=0, head_index=1, stat="SCD")
    df = df.rename(columns={"SCD_m0": "SCD_m0_bg0"})
    df["INS-16_m0_bg0"] = average_stat_over_targets(
        df, model_index=0, head_index=1, stat="INS-16"
    )["INS-16_m0"]
    df["INS-64_m0_bg0"] = average_stat_over_targets(
        df, model_index=0, head_index=1, stat="INS-64"
    )["INS-64_m0"]
    df = df.drop(
        columns=[f"SCD_h1_m0_t{x}" for x in range(6)]
        + [f"INS-16_h1_m0_t{x}" for x in range(6)]
        + [f"INS-64_h1_m0_t{x}" for x in range(6)]
        + ["background_index"]
    )

    # Process and add data for all models and backgrounds
    for model_index in range(4):
        start_from = 1 if model_index == 0 else 0
        for bg_index in range(start_from, 10):
            stats = process_data(model_index, bg_index)
            for stat, df_stat in stats.items():
                df[f"{stat}_m{model_index}_bg{bg_index}"] = df_stat[
                    f"{stat}_m{model_index}_bg{bg_index}"
                ]

    # Averaging
    for stat in ["SCD", "INS-16", "INS-64"]:
        df[stat] = df[
            [
                f"{stat}_m{model_ind}_bg{bg_ind}"
                for model_ind in range(4)
                for bg_ind in range(10)
            ]
        ].mean(axis=1)

    return df
