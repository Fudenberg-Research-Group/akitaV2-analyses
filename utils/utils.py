import numpy as np
import pandas as pd
import pysam
import scipy

from akita_utils.format_io import h5_to_df
from akita_utils.dna_utils import (dna_rc, dna_1hot_index, dna_1hot)

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
        
    df[f"{stat}_m{model_index}"] = df[[f"{stat}_h{head_index}_m{model_index}_t{target_index}" for target_index in range(target_indices)]].mean(axis=1)
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


def average_stat_over_backgrounds(df,
                                model_index=0,
                                head_index=1,
                                num_backgrounds=10,
                                stat="SCD",
                                columns_to_keep = ["chrom", "end", "start", "strand", "seq_id"],
                                keep_background_columns=True):
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

    num_sites = len(df) // 10
    output_df = df[columns_to_keep][:num_sites]

    for bg_index in range(num_backgrounds):
        output_df[f"{stat}_bg{bg_index}"] = df[df["background_index"] == bg_index][f"{stat}_m{model_index}"].values

    output_df[f"{stat}_m{model_index}"] = output_df[[f"{stat}_bg{bg_index}" for bg_index in range(num_backgrounds)]].mean(axis=1)

    if keep_background_columns == False:
        output_df = output_df.drop(columns=[f"{stat}_bg{bg_index}" for bg_index in range(num_backgrounds)])
    
    return output_df


def read_and_average_virtual_exp(data_dir, stat_to_average="SCD", all_calculated_stats=["SCD", "INS-16", "INS-64"], model_numbers=8):
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
    df_m0 = h5_to_df(data_dir+"/model_0.h5", all_calculated_stats, average=False) 
    df_m1 = h5_to_df(data_dir+"/model_1.h5", all_calculated_stats, average=False) 
    df_m2 = h5_to_df(data_dir+"/model_2.h5", all_calculated_stats, average=False) 
    df_m3 = h5_to_df(data_dir+"/model_3.h5", all_calculated_stats, average=False) 
    df_m4 = h5_to_df(data_dir+"/model_4.h5", all_calculated_stats, average=False) 
    df_m5 = h5_to_df(data_dir+"/model_5.h5", all_calculated_stats, average=False) 
    df_m6 = h5_to_df(data_dir+"/model_6.h5", all_calculated_stats, average=False) 
    df_m7 = h5_to_df(data_dir+"/model_7.h5", all_calculated_stats, average=False) 

    print("averaging over targets")
    df_m0_tg = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat=stat_to_average)
    df_m1_tg = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat=stat_to_average)
    df_m2_tg = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat=stat_to_average)
    df_m3_tg = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat=stat_to_average)
    df_m4_tg = average_stat_over_targets(df_m4, model_index=4, head_index=1, stat=stat_to_average)
    df_m5_tg = average_stat_over_targets(df_m5, model_index=5, head_index=1, stat=stat_to_average)
    df_m6_tg = average_stat_over_targets(df_m6, model_index=6, head_index=1, stat=stat_to_average)
    df_m7_tg = average_stat_over_targets(df_m7, model_index=7, head_index=1, stat=stat_to_average)

    print("averaging over backgrounds")
    df_m0_tgbg = average_stat_over_backgrounds(df_m0_tg, model_index=0, head_index=1, stat=stat_to_average)
    df_m1_tgbg = average_stat_over_backgrounds(df_m1_tg, model_index=1, head_index=1, stat=stat_to_average)
    df_m2_tgbg = average_stat_over_backgrounds(df_m2_tg, model_index=2, head_index=1, stat=stat_to_average)
    df_m3_tgbg = average_stat_over_backgrounds(df_m3_tg, model_index=3, head_index=1, stat=stat_to_average)
    df_m4_tgbg = average_stat_over_backgrounds(df_m4_tg, model_index=4, head_index=1, stat=stat_to_average)
    df_m5_tgbg = average_stat_over_backgrounds(df_m5_tg, model_index=5, head_index=1, stat=stat_to_average)
    df_m6_tgbg = average_stat_over_backgrounds(df_m6_tg, model_index=6, head_index=1, stat=stat_to_average)
    df_m7_tgbg = average_stat_over_backgrounds(df_m7_tg, model_index=7, head_index=1, stat=stat_to_average)
    
    print(f"collecting data for {stat_to_average}")
    stat_collected = pd.concat([df_m0_tgbg["chrom"], df_m0_tgbg["end"], df_m0_tgbg["start"], df_m0_tgbg["strand"], df_m0_tgbg[f"{stat_to_average}_m0"], df_m1_tgbg[f"{stat_to_average}_m1"], df_m2_tgbg[f"{stat_to_average}_m2"], df_m3_tgbg[f"{stat_to_average}_m3"], df_m4_tgbg[f"{stat_to_average}_m4"], df_m5_tgbg[f"{stat_to_average}_m5"], df_m6_tgbg[f"{stat_to_average}_m6"], df_m7_tgbg[f"{stat_to_average}_m7"]], axis=1)

    # averaging over models
    stat_collected[f"{stat_to_average}"] = stat_collected[[f"{stat_to_average}_m{model_index}" for model_index in range(model_numbers)]].mean(axis=1)

    return stat_collected


def read_and_average_genomic_exp(data_dir, stat_to_average="SCD", all_calculated_stats=["SCD", "INS-16", "INS-64"], model_numbers=8):
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
    df_m0 = h5_to_df(data_dir+"/model_0.h5", all_calculated_stats, average=False) 
    df_m1 = h5_to_df(data_dir+"/model_1.h5", all_calculated_stats, average=False) 
    df_m2 = h5_to_df(data_dir+"/model_2.h5", all_calculated_stats, average=False) 
    df_m3 = h5_to_df(data_dir+"/model_3.h5", all_calculated_stats, average=False) 
    df_m4 = h5_to_df(data_dir+"/model_4.h5", all_calculated_stats, average=False) 
    df_m5 = h5_to_df(data_dir+"/model_5.h5", all_calculated_stats, average=False) 
    df_m6 = h5_to_df(data_dir+"/model_6.h5", all_calculated_stats, average=False) 
    df_m7 = h5_to_df(data_dir+"/model_7.h5", all_calculated_stats, average=False) 

    print("averaging over targets")
    df_m0_tg = average_stat_over_targets(df_m0, model_index=0, head_index=1, stat=stat_to_average)
    df_m1_tg = average_stat_over_targets(df_m1, model_index=1, head_index=1, stat=stat_to_average)
    df_m2_tg = average_stat_over_targets(df_m2, model_index=2, head_index=1, stat=stat_to_average)
    df_m3_tg = average_stat_over_targets(df_m3, model_index=3, head_index=1, stat=stat_to_average)
    df_m4_tg = average_stat_over_targets(df_m4, model_index=4, head_index=1, stat=stat_to_average)
    df_m5_tg = average_stat_over_targets(df_m5, model_index=5, head_index=1, stat=stat_to_average)
    df_m6_tg = average_stat_over_targets(df_m6, model_index=6, head_index=1, stat=stat_to_average)
    df_m7_tg = average_stat_over_targets(df_m7, model_index=7, head_index=1, stat=stat_to_average)
    
    print(f"collecting data for {stat_to_average}")
    stat_collected = pd.concat([df_m0_tg["chrom"], df_m0_tg["end"], df_m0_tg["start"], df_m0_tg["strand"], df_m0_tg[f"{stat_to_average}_m0"], df_m1_tg[f"{stat_to_average}_m1"], df_m2_tg[f"{stat_to_average}_m2"], df_m3_tg[f"{stat_to_average}_m3"], df_m4_tg[f"{stat_to_average}_m4"], df_m5_tg[f"{stat_to_average}_m5"], df_m6_tg[f"{stat_to_average}_m6"], df_m7_tg[f"{stat_to_average}_m7"]], axis=1)

    # averaging over models
    stat_collected[f"{stat_to_average}"] = stat_collected[[f"{stat_to_average}_m{model_index}" for model_index in range(model_numbers)]].mean(axis=1)

    return stat_collected

def collect_flanked_sequences(sites, flank_length=30, genome_path="/project/fudenber_735/genomes/mm10/mm10.fa"):
    """
    Extracts and processes DNA sequences from a given genome, flanking specified genomic sites.

    This function reads a genome from a specified path and retrieves sequences that flank genomic sites of interest. 
    Each site is extended by a given flank length on both sides. The function also accounts for the strand 
    orientation, providing the reverse complement sequence if the strand is negative.

    Parameters:
    - sites (DataFrame): A pandas DataFrame containing the genomic sites of interest. The DataFrame must 
      include columns 'chrom', 'start', 'end', and 'strand'.
    - flank_length (int, optional): The number of base pairs to include on each side of the site. Default is 30.
    - genome_path (str, optional): The file path to the genome fasta file. Default is "/project/fudenber_735/genomes/mm10/mm10.fa".

    Returns:
    - numpy.ndarray: An array where each element is a one-hot encoded representation of the flanked sequence at each site.
    """
    genome_open = pysam.Fastafile(genome_path)
    sites_dna_num = []
    
    for i in range(len(sites)):
        chrm, start, end, strand = sites.iloc[i][["chrom","start","end","strand"]]
        start = start-flank_length
        end = end+flank_length
        seq =  genome_open.fetch(chrm, start, end).upper()
        if strand == '-':
            seq = dna_rc(seq)
        sites_dna_num.append(dna_1hot_index(seq))
    
    genome_open.close()
    sites_dna_num = np.array(sites_dna_num)
    
    return sites_dna_num

def reorder_by_hamming_dist(dna_matrix, sub_index=(0,-1)):    
    """
    Reorders a matrix of DNA sequences based on their pairwise Hamming distances.

    This function calculates the pairwise Hamming distances between rows (DNA sequences) in a given matrix. 
    It then reorders the matrix such that similar sequences (based on the Hamming distance) are grouped together. 
    The function allows for considering only a subset of each sequence for the distance calculation.

    Parameters:
    - dna_matrix (numpy.ndarray): A 2D numpy array where each row represents a DNA sequence.
    - sub_index (tuple, optional): A tuple (start, end) indicating the subset of each sequence to consider 
      for the Hamming distance calculation. Default is (0, -1) which considers the full length of each sequence.

    Returns:
    - numpy.ndarray: The reordered matrix of DNA sequences, where sequences are ordered based on their 
      similarity (lower Hamming distance).

    Note: This function assumes that all sequences in dna_matrix are of equal length.
    """
    num_seqs = len(dna_matrix)   
    seq_dist = np.zeros((num_seqs,num_seqs))
    
    for i in range(num_seqs):
        seq_i = dna_matrix[i][sub_index[0]:sub_index[1]]
        for j in range(num_seqs):
            if i<j:
                seq_j = dna_matrix[j][sub_index[0]:sub_index[1]]
                seq_dist[i,j] = 1 - scipy.spatial.distance.hamming(seq_i,seq_j)
    seq_dist = seq_dist + seq_dist.T           
    
    reording = scipy.cluster.hierarchy.leaves_list(
                    scipy.cluster.hierarchy.linkage(seq_dist))
    return dna_matrix[reording]


def prepare_nt_count_table(sites, flank_length=30, genome_path="/project/fudenber_735/genomes/mm10/mm10.fa"):
    """
    Prepares a nucleotide count table for a set of genomic sites with specified flanking regions.

    This function processes a list of genomic sites and extracts the corresponding DNA sequences, 
    including flanking regions, from a specified genome. It then counts the occurrences of each 
    nucleotide at each position across all sequences and compiles these counts into a table.

    Parameters:
    - sites (DataFrame): A pandas DataFrame containing the genomic sites of interest. The DataFrame must 
      include columns 'chrom', 'start', 'end', and 'strand'.
    - flank_length (int, optional): The number of base pairs to include on each side of the site. Default is 30.
    - genome_path (str, optional): The file path to the genome fasta file. Default is "/project/fudenber_735/genomes/mm10/mm10.fa".

    Returns:
    - numpy.ndarray: A 2D numpy array with shape (sequence_length, 4), where sequence_length is the length 
      of the flanked site. Each column corresponds to one of the four nucleotides (A, C, G, T), and each 
      row represents the count of each nucleotide at that position across all sequences.
    """
    seq_length = flank_length * 2 + (sites.iloc[0]["end"] - sites.iloc[0]["start"])
    genome_open = pysam.Fastafile(genome_path)
    
    nt_count = np.zeros(shape=(seq_length, 4))
    for i in range(len(sites)):
        chrm, start, end, strand = sites.iloc[i][["chrom","start","end","strand"]]
        start = start-flank_length
        end = end+flank_length
        seq =  genome_open.fetch(chrm, start, end).upper()
        if strand == '-':
            seq = dna_rc(seq)
        nt_count = nt_count + dna_1hot(seq)
    genome_open.close()
    return nt_count
