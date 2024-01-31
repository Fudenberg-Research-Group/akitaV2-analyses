import numpy as np
import pandas as pd
import pysam
import scipy

from akita_utils.dna_utils import dna_rc, dna_1hot_index, dna_1hot
from akita_utils.tsv_gen_utils import filter_dataframe_by_column

########################################
# flanking sequences analysis          #
########################################

def collect_flanked_sequences(
    sites,
    flank_length=30,
    genome_path="/project/fudenber_735/genomes/mm10/mm10.fa",
):
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
        chrm, start, end, strand = sites.iloc[i][
            ["chrom", "start", "end", "strand"]
        ]
        start = start - flank_length
        end = end + flank_length
        seq = genome_open.fetch(chrm, start, end).upper()
        if strand == "-":
            seq = dna_rc(seq)
        sites_dna_num.append(dna_1hot_index(seq))

    genome_open.close()
    sites_dna_num = np.array(sites_dna_num)

    return sites_dna_num


def reorder_by_hamming_dist(dna_matrix, sub_index=(0, -1)):
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
    seq_dist = np.zeros((num_seqs, num_seqs))

    for i in range(num_seqs):
        seq_i = dna_matrix[i][sub_index[0] : sub_index[1]]
        for j in range(num_seqs):
            if i < j:
                seq_j = dna_matrix[j][sub_index[0] : sub_index[1]]
                seq_dist[i, j] = 1 - scipy.spatial.distance.hamming(
                    seq_i, seq_j
                )
    seq_dist = seq_dist + seq_dist.T

    reording = scipy.cluster.hierarchy.leaves_list(
        scipy.cluster.hierarchy.linkage(seq_dist)
    )
    return dna_matrix[reording]


def prepare_nt_count_table(
    sites,
    flank_length=30,
    genome_path="/project/fudenber_735/genomes/mm10/mm10.fa",
):
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
    seq_length = flank_length * 2 + (
        sites.iloc[0]["end"] - sites.iloc[0]["start"]
    )
    genome_open = pysam.Fastafile(genome_path)

    nt_count = np.zeros(shape=(seq_length, 4))
    for i in range(len(sites)):
        chrm, start, end, strand = sites.iloc[i][
            ["chrom", "start", "end", "strand"]
        ]
        start = start - flank_length
        end = end + flank_length
        seq = genome_open.fetch(chrm, start, end).upper()
        if strand == "-":
            seq = dna_rc(seq)
        nt_count = nt_count + dna_1hot(seq)
    genome_open.close()
    return nt_count


########################################
# DataFrame functions                  #
########################################


def split_by_percentile_groups(df, column_to_split, num_classes, 
                               upper_percentile=100, lower_percentile=0, 
                               category_colname="category"):
    """
    Splits a dataframe into distinct groups based on the percentile ranges of a specified column.
    Each group represents a percentile range based on the number of classes specified.

    Parameters
    ----------
    df : DataFrame
        The input pandas dataframe.
    column_to_split : str
        The column based on which the dataframe is split into percentile groups.
    num_classes : int
        The number of classes to split the dataframe into.
    upper_percentile : int, default 100
        The upper limit of the percentile range. Typically set to 100.
    lower_percentile : int, default 0
        The lower limit of the percentile range. Typically set to 0.
    category_colname : str, default "category"
        The name of the new column to be added to the dataframe, indicating the category of each row based on percentile range.

    Returns
    -------
    DataFrame
        A new dataframe with an additional column named as specified by 'category_colname'.
        This column contains categorical labels corresponding to the specified percentile ranges.
    """
    bounds = np.linspace(lower_percentile, (upper_percentile-lower_percentile), num_classes + 1, dtype="int")
    df_out = pd.DataFrame()

    for i in range(num_classes):
        
        group_df = filter_dataframe_by_column(
            df,
            column_name=column_to_split,
            upper_threshold=bounds[i+1],
            lower_threshold=bounds[i],
            drop_duplicates=False
        )
        group_df[category_colname] = f"Group_{i}"
        df_out = pd.concat([df_out, group_df])

    return df_out


