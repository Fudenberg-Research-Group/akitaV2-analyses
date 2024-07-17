# This script prepares a dataframe specifying genomic sequence loci for shuffling experiments.
# It reads genomic intervals from a BED file and calculates GC content using a reference genome fasta file.
# Parameters such as shuffle sizes, motif detection thresholds, mutation methods, output filename,
# sample size, and locus selection mode are configurable via command-line arguments.
# The script generates all possible combinations of these parameters and saves them as a TSV file.
#
# Inputs:
# - Genome fasta file containing the reference genome sequences.
# - BED file specifying genomic intervals to prepare for shuffling.
#
# Parameters:
# - shuffle_parameter: K-mer sizes for sequence shuffling.
# - ctcf_detection_threshold: Accuracy thresholds for motif identification.
# - mutation_method: Methods for sequence mutation (e.g., permutation, randomization, motif masking).
# - output_filename: Name of the output TSV file.
# - num_seqs: Number of sequences to sample from the BED file based on GC content.
# - mode: Criteria for selecting genomic loci based on GC content (uniform, tail, head, random).
#
# Output:
# - A TSV file with all parameter combinations and associated sequence loci information.
#
# Usage notes:
# - Multiple values for each parameter can be specified via the command-line interface (e.g., --shuffle_parameter 2 4 8).
# - Example BED files for mouse (mm10) and human (hg38) genomes are provided as reference.
#   Adjust paths accordingly for different datasets.
#
# Example command-line usage:
# python generate_shuffled_seqs_df.py -f genome.fasta -seq_bed_file intervals.bed --shuffle_parameter 2 4 8 --output_filename shuffled_seqs.tsv --num_seqs 20 --mode uniform

# import general libraries
import itertools
import pandas as pd
import bioframe
import argparse
from akita_utils.tsv_utils import filter_dataframe_by_column


################################################################################
# main
################################################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="genome_fasta", help="fasta file", required=True)
    parser.add_argument(
        "-seq_bed_file",
        dest="seq_bed_file",
        help="bed file for the seqs under investigation",
        required=True,
    )
    parser.add_argument(
        "--output_filename",
        dest="output_filename",
        default="data/shuffled_seqs.tsv",
        help="output_filename",
    )
    parser.add_argument(
        "--shuffle_parameter",
        nargs="+",
        default=[2, 4, 8],
        type=int,
        help="List of integers sepaerated by spaces eg 2 4",
    )
    parser.add_argument(
        "--ctcf_detection_threshold",
        default=[8],
        nargs="+",
        type=int,
        help="threshold of (CTCF PWM) * DNA OHE window value",
    )
    parser.add_argument(
        "--mutation_method",
        nargs="+",
        default=["permute_whole_seq"],
        help="mutation methods from ['permute_whole_seq','randomise_whole_seq','randomise_motif','permute_motif','mask_motif'], space seperated",
    )
    parser.add_argument(
        "--num_seqs",
        type=int,
        default=20,
        help="number of seqs to select from dataframe",
    )
    parser.add_argument("--mode", default="uniform", help="loci selection criteria")
    args = parser.parse_args()

    # prepare dataframe with chromosomes and calculate GC content(using bioframe)
    seq_df = pd.read_csv(
        args.seq_bed_file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "fold"],
    )
    general_seq_gc_df = bioframe.frac_gc(seq_df, bioframe.load_fasta(args.genome_fasta), return_input=True)

    grid_search_params = {
        "shuffle_parameter": args.shuffle_parameter,
        "ctcf_detection_threshold": args.ctcf_detection_threshold,
        "mutation_method": args.mutation_method,
    }
    
    # sampling seq_df dataframe respecting GC content
    seq_gc_df = filter_dataframe_by_column(
        general_seq_gc_df,
        column_name="GC",
        upper_threshold=99,
        lower_threshold=1,
        filter_mode=args.mode,
        num_rows=args.num_seqs,
    )

    # fixing locus specific chacteristics together before grid_search
    seq_gc_df = (
        seq_gc_df["chrom"].map(str)
        + ","
        + seq_gc_df["start"].map(str)
        + ","
        + seq_gc_df["end"].map(str)
        + ","
        + seq_gc_df["GC"].map(str)
    )
    locus_list = seq_gc_df.values.tolist()

    grid_search_params["locus_specification"] = locus_list

    grid_search_param_set = list(itertools.product(*[v for v in grid_search_params.values()]))
    parameters_combo_dataframe = pd.DataFrame(grid_search_param_set, columns=grid_search_params.keys())
    parameters_combo_dataframe[["chrom", "start", "end", "GC_content"]] = parameters_combo_dataframe[
        "locus_specification"
    ].str.split(",", expand=True)
    parameters_combo_dataframe = parameters_combo_dataframe.drop(columns=["locus_specification"])
    
    parameters_combo_dataframe.to_csv(f"{args.output_filename}", sep="\t", index=False)


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
    