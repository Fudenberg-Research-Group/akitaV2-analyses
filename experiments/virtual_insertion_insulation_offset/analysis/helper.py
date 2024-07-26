import numpy as np
import matplotlib.pyplot as plt

plt.style.use(
    "/home1/smaruj/akitaX1-analyses/figures/plot_styles/global_plotting_style.mplstyle"
)

from akita_utils.format_io import h5_to_df
from akita_utils.df_utils import average_stat_over_targets


def read_and_average_insulation_offset_data(
    data_dir, keys_to_ignore=["insertion_SCD", "disruption_SCD"]
):
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
    stats = ["OFF-16", "OFF-64", "OFF-128"]  # Statistics to be averaged

    # Read data for each model into a dictionary of DataFrames
    dfs = {
        model: h5_to_df(
            f"{data_dir}model_{model}.h5",
            ["SCD"] + stats,
            ignore_keys=keys_to_ignore,
            average=False,
        )
        for model in models
    }

    # Initialize DataFrame to store results
    # df = pd.DataFrame()
    df = dfs[0][["chrom", "start", "end", "strand", "orientation"]].copy()

    # Calculate averages for each statistic across all models
    for stat in stats:
        for model in models:
            df[f"{stat}_m{model}"] = average_stat_over_targets(
                dfs[model], model_index=model, head_index=1, stat=stat
            )[f"{stat}_m{model}"]

        # Compute overall average for each statistic
        df[stat] = df[[f"{stat}_m{model}" for model in models]].mean(axis=1)

    return df


def extract_window_from_vector(vector, window=10, width=3):
    """
    Extract a window of size 3*window centered around the center of the vector.

    Parameters
    ------------
    vector : list or numpy array
        Input vector.
    window : int, optional
        Size of the window. Default is 16.

    Returns
    ---------
    window_vector : list or numpy array
        Extracted window of size 3*window.
    """
    center_index = len(vector) // 2
    start_index = max(center_index - width * window, 0)
    end_index = min(center_index + width * window, len(vector))
    window_vector = vector[start_index:end_index]
    return window_vector


def slide_diagonal_insulation(target_map, window=10):
    """
    Calculate insulation by sliding a window along the diagonal of the target map.

    Parameters
    ------------
    target_map : numpy array
        Array with contact change maps predicted by Akita, usually of a size (512 x 512)
    window : int
        Size of the sliding window.

    Returns
    ---------
    scores : list
        List of ISN-window scores for each position along the diagonal.
    """

    map_size = target_map.shape[0]
    scores = np.empty((map_size,))
    scores[:] = np.nan

    for mid in range(window, (map_size - window)):
        lo = max(0, mid + 1 - window)
        hi = min(map_size, mid + window)
        score = np.nanmean(target_map[lo : (mid + 1), mid:hi])
        scores[mid] = score

    return scores


def min_insulation_offset_from_center(
    target_map, window=10, crop_around_center=True, crop_width=3
):
    """
    Calculate the offset from the center position for the position along the diagonal
    with the minimum insulation score for a sliding window along the diagonal.

    Parameters
    ------------
    target_map : numpy array
        Array with contact change maps predicted by Akita, usually of a size (512 x 512).
    window : int, optional
        Size of the sliding window for insulation calculation. Default is 10.

    Returns
    ---------
    offset_from_center : int
        Offset from the center position based on the position with the minimum insulation score.
    """
    map_size = target_map.shape[0]
    center_position = map_size // 2
    bin_shift = 0

    insulation_scores = slide_diagonal_insulation(target_map, window)

    if crop_around_center:
        bin_shift = max(center_position - crop_width * window, 0)
        insulation_scores = extract_window_from_vector(
            insulation_scores, window=window, width=3
        )

    min_score = np.nanmin(insulation_scores)

    # Find indices with min_score
    indices_tuple = np.where(insulation_scores == min_score)
    indices_list = list(indices_tuple[0])

    # in case there are more than one index with min_score
    # Calculate midpoint of the array
    midpoint = len(insulation_scores) // 2

    # Find index closest to the midpoint among indices with the same score
    closest_index = None
    min_difference = float("inf")
    for index in indices_list:
        difference = abs(index - midpoint)
        if difference < min_difference:
            closest_index = index
            min_difference = difference

    closest_index = closest_index + bin_shift
    offset_from_center = closest_index - center_position
    return offset_from_center


def calculate_offset_INS(
    map_matrix, window=10, crop_around_center=True, crop_width=3
):
    """
    Calculate insulation in a window-size diamond around the central pixel
    for a set of num_targets contact difference maps.

    Parameters
    ------------
    map_matrix : numpy array
        Array with contact change maps predicted by Akita, usually of a size (512 x 512 x num_targets)
    window : int
        Size of a diamond taken into account in a metric calculation.

    Returns
    ---------
    scores : num_targets-long vector with ISN-window offset
    """

    num_targets = map_matrix.shape[-1]
    scores = np.zeros((num_targets,))
    for target_index in range(num_targets):
        offset_from_center = min_insulation_offset_from_center(
            map_matrix[:, :, target_index],
            window=window,
            crop_around_center=crop_around_center,
            crop_width=crop_width,
        )
        scores[target_index] = offset_from_center
    return scores


def plot_insulation_scores(
    insulation_scores_ins16,
    insulation_scores_ins64,
    insulation_scores_ins128,
    window_to_plot=40,
    save_path=None,
):
    """
    Plot insulation scores against positions along the diagonal.

    Parameters
    ------------
    insulation_scores : list
        List of insulation scores.

    Returns
    ---------
    None
    """
    midpoint = len(insulation_scores_ins16) // 2
    positions_ins16 = [i for i in range(2 * window_to_plot)]

    # changing windows from bins to kbs
    window_128 = 128 * 2048 // 1000
    window_64 = 64 * 2048 // 1000
    window_16 = 16 * 2048 // 1000

    plt.figure(figsize=(8, 5))
    plt.plot(
        positions_ins16,
        insulation_scores_ins128[
            midpoint - window_to_plot : midpoint + window_to_plot
        ],
        marker="o",
        linestyle="-",
        c="#7570b3",
        label=f"{window_128} kb",
    )
    plt.plot(
        positions_ins16,
        insulation_scores_ins64[
            midpoint - window_to_plot : midpoint + window_to_plot
        ],
        marker="o",
        linestyle="-",
        c="#d95f02",
        label=f"{window_64} kb",
    )
    plt.plot(
        positions_ins16,
        insulation_scores_ins16[
            midpoint - window_to_plot : midpoint + window_to_plot
        ],
        marker="o",
        linestyle="-",
        c="#1b9e77",
        label=f"{window_16} kb",
    )

    # Set global y-axis range starting from -0.7
    plt.ylim(-0.75, 0)
    # Adjust the upper limit based on maximum score

    plt.xlabel("Middle position along diagonal")
    plt.ylabel("Insulation Score")

    # Custom x-ticks
    tick_positions = range(
        0, 2 * window_to_plot + 10, 10
    )  # Positions for ticks every 10 units
    tick_labels = [
        str(i + midpoint - window_to_plot) for i in tick_positions
    ]  # Adjusted labels starting from midpoint-window_to_plot

    plt.xticks(tick_positions, tick_labels)  # Apply custom ticks and labels

    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)

    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    else:
        plt.show()
