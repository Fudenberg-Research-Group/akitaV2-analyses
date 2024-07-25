import pandas as pd
from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import average_stat_for_shift


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