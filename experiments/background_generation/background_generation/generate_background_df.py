# import general libraries

import itertools
import pandas as pd
import bioframe
import argparse
from akita_utils.tsv_utils import filter_dataframe_by_column


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        dest="genome_fasta",
        help="fasta file",
        required=True,
    )
    parser.add_argument(
        "-seq_bed_file",
        dest="seq_bed_file",
        help="bed file for the seqs under investigation",
        required=True,
    )
    parser.add_argument(
        "--output_filename",
        dest="output_filename",
        default="data/flat_seqs_test.tsv",
        help="output_filename",
    )
    parser.add_argument(
        "--shuffle_parameter",
        nargs="+",
        default=[8],
        type=int,
        help="integers separated by spaces",
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
        "--num_backgrounds",
        type=int,
        default=10,
        help="number of loci to select for background creation",
    )
    parser.add_argument(
        "--mode",
        default="uniform",
        help="loci selection criteria",
    )
    parser.add_argument(
        "--SCD_threshold",
        type=float,
        nargs="+",
        default=[40],
        help="maximum allowable map score, SCD",
    )

    args = parser.parse_args()

    # prepare dataframe with chromosomes and calculate GC content(using bioframe)
    seq_df = pd.read_csv(
        args.seq_bed_file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "fold"],
    )
    general_seq_gc_df = bioframe.frac_gc(
        seq_df, bioframe.load_fasta(args.genome_fasta), return_input=True
    )

    grid_params = {
        "shuffle_parameter": args.shuffle_parameter,
        "ctcf_detection_threshold": args.ctcf_detection_threshold,
        "mutation_method": args.mutation_method,
        "map_score_threshold": args.SCD_threshold,
    }

    # sampling seq_df dataframe respecting GC content
    seq_gc_df = filter_dataframe_by_column(
        general_seq_gc_df,
        column_name="GC",
        upper_threshold=99,
        lower_threshold=1,
        filter_mode=args.mode,
        num_rows=args.num_backgrounds,
    )

    # fixing locus specific variables together before grid creation
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

    grid_params["locus_specification"] = locus_list

    grid_param_set = list(
        itertools.product(*[v for v in grid_params.values()])
    )
    parameters_combo_dataframe = pd.DataFrame(
        grid_param_set, columns=grid_params.keys()
    )
    parameters_combo_dataframe[
        ["chrom", "start", "end", "GC_content"]
    ] = parameters_combo_dataframe["locus_specification"].str.split(
        ",", expand=True
    )
    parameters_combo_dataframe = parameters_combo_dataframe.drop(
        columns=["locus_specification"]
    )

    parameters_combo_dataframe.to_csv(
        f"{args.output_filename}", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
