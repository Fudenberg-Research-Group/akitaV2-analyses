import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from collections import Counter
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


def load_sequences(fasta_file, start=0, end=-1):
    sequences = []
    with open(fasta_file, "r") as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence[start:end])
                    sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence[start:end])
    return sequences

def load_top_bottom_sequences(fasta_path):
    whole = load_sequences(fasta_path, start=0, end=-1)
    upstream = load_sequences(fasta_path, start=0, end=30)
    downstream = load_sequences(fasta_path, start=-30, end=None)
    return whole, upstream, downstream


def calculate_gc_content_at_positions(sequences):
    # Find the length of the shortest sequence
    min_length = min(len(seq) for seq in sequences)

    # Initialize counts arrays
    gc_counts = np.zeros(min_length)
    total_counts = np.zeros(min_length)

    for sequence in sequences:
        # Trim sequence to the length of the shortest sequence
        trimmed_seq = sequence[:min_length]
        for i, nucleotide in enumerate(trimmed_seq):
            if nucleotide in "GC":
                gc_counts[i] += 1
            total_counts[i] += 1

    gc_content_at_positions = (gc_counts / total_counts) * 100
    return gc_content_at_positions


def count_kmers_total(sequences, k):
    total_kmer_counts = Counter()
    for sequence in sequences:
        k_mers = generate_kmers(sequence, k)
        total_kmer_counts.update(k_mers)
    return total_kmer_counts


def generate_kmers(sequence, k):
    k_mers = []
    for i in range(len(sequence) - k + 1):
        k_mers.append(sequence[i : i + k])
    return k_mers


def count_kmers(sequences, k):
    kmer_counts = []
    for sequence in sequences:
        k_mers = generate_kmers(sequence, k)
        kmer_counts.append(Counter(k_mers))
    return kmer_counts


def create_kmer_vector(kmer_counts, k):
    all_kmers = sorted(set(kmer for counts in kmer_counts for kmer in counts))
    vectors = []
    for counts in kmer_counts:
        vector = [counts[kmer] for kmer in all_kmers]
        vectors.append(vector)
    return np.array(vectors), all_kmers


def apply_pca(kmer_vectors, n_components=2):
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(kmer_vectors)
    return pca_result


def process_sequences(sequences, k):
    counts = {key: count_kmers(seq, k=k) for key, seq in sequences.items()}
    vectors = {key: create_kmer_vector(counts[key], k=k)[0] for key in counts}
    pca_results = {key: apply_pca(vectors[key]) for key in vectors}
    return pca_results


def read_meme_pwm_as_numpy(filename):
    pwm_list = []  # List to store PWM rows
    
    with open(filename, 'r') as file:
        in_matrix_section = False
        
        for line in file:
            line = line.strip()
            
            # Check if we are reading the PWM matrix
            if line.startswith("letter-probability matrix"):
                in_matrix_section = True  # Start reading matrix data
                continue  # Skip this header line
            
            # If we are in the matrix section, process the rows
            if in_matrix_section and line:
                pwm_row = [float(value) for value in line.split()]  # Parse values
                pwm_list.append(pwm_row)  # Append to the PWM list
            
            # If we encounter a new MOTIF or the end of file, stop matrix reading
            if line.startswith("MOTIF") and in_matrix_section:
                break
    
    # Convert the list to a numpy array
    pwm_array = np.array(pwm_list)
    
    return pwm_array


# Weighted difference of probabilities
# based on:
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0767-x

def jensen_shannon_divergence(p, q):
    """
    Compute the Jensen-Shannon divergence between two probability distributions p and q
    """
    m = 0.5 * (p + q)
    divergence = 0.5 * np.sum(p * np.log2(np.divide(p, m, out=np.zeros_like(p), where=m != 0)), axis=1) + \
                 0.5 * np.sum(q * np.log2(np.divide(q, m, out=np.zeros_like(q), where=m != 0)), axis=1)
    return divergence

def calculate_weight(p, q):
    """
    Compute the weight r_ell,a as described in the method.
    """
    abs_diff = np.abs(p - q)
    sum_abs_diff = np.sum(abs_diff, axis=1, keepdims=True)
    
    # To avoid division by zero, use np.where to handle p == q
    weight = np.divide(p - q, sum_abs_diff, out=np.zeros_like(p), where=sum_abs_diff != 0)
    
    return weight

def weighted_difference_of_probabilities(pwm1, pwm2):
    """
    Compute the weighted difference of probabilities using Jensen-Shannon divergence
    between two PWM matrices of size (n, 4).
    
    Parameters:
    pwm1: First PWM matrix of shape (n, 4)
    pwm2: Second PWM matrix of shape (n, 4)
    
    Returns:
    H_ell_a: Weighted heights for each symbol a at each position ell
    """
    # Step 1: Calculate Jensen-Shannon divergence H_ell
    H_ell = jensen_shannon_divergence(pwm1, pwm2)
    
    # Step 2: Calculate weights r_ell,a
    weights = calculate_weight(pwm1, pwm2)
    
    # Step 3: Calculate weighted heights H_ell,a
    H_ell_a = weights * H_ell[:, np.newaxis]  # Broadcasting H_ell to match the shape of weights
    
    H_ell_a = np.nan_to_num(H_ell_a, nan=0.0)
    
    return H_ell_a


def normalize_pwm(pwm, epsilon=1e-10):
    """
    Add a small value (epsilon) to the PWM matrix to avoid zeros, then renormalize
    so that each row sums to 1.
    
    Parameters:
    pwm: Input PWM matrix of shape (n, 4)
    epsilon: Small value to add to each element in the PWM matrix
    
    Returns:
    normalized_pwm: Renormalized PWM matrix of shape (n, 4)
    """
    pwm_with_epsilon = pwm + epsilon
    row_sums = np.sum(pwm_with_epsilon, axis=1, keepdims=True)
    normalized_pwm = pwm_with_epsilon / row_sums  # Renormalize each row to sum to 1
    
    return normalized_pwm