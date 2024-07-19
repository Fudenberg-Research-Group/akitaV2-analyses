# Description:
# This script processes a TSV file containing CTCF-binding site coordinates by generating orientations, assigning background indices, 
# and adding flank and spacer information. The processed data is then saved to a new TSV file.
#
# Inputs:
# - --input-tsv-file: Path to the input TSV file with CTCF-binding site coordinates [Default: /home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv].
# - --num-inserts: Number of inserts (core + flanks) to be inserted into a background sequence [Default: 6].
# - --flank-length: Range of right and left flank to be tested [Default: 30].
# - --flank-spacer-sum: Sum of flank and spacer distances to keep distances between CTCFs constant [Default: 90].
# - --backgrounds-indices: Comma-separated indices of background sequences [Default: "0,1,2,3,4,5,6,7,8,9"].
# - --num-ctcf-sites: Number of CTCF sites to be tested [Default: 100].
# - --output-filename: Filename for the output TSV file [Default: "out.tsv"].
#
# Outputs:
# - TSV file with additional columns for orientation, background index, and flank/spacer information.
#
# Example command-line usage:
# python generate_insulation_oddset_orientation_df.py --input-tsv-file input.tsv --num-inserts 6 --flank-length 30 --flank-spacer-sum 90 --backgrounds-indices "0,1,2,3" --num-ctcf-sites 100 --output-filename processed_ctcfs.tsv

from optparse import OptionParser
import pandas as pd
import numpy as np

from cooltools.lib import numutils

from akita_utils.tsv_utils import (
    add_orientation,
    add_background,
    add_diff_flanks_and_const_spacer,
)

################################################################################
# helper function
################################################################################

def generate_orientations(num_sites):
    orientations = []
    
    # Generate symmetric orientations
    symmetric_orientation = ">" * num_sites
    orientations.append(symmetric_orientation)
    
    # Generate asymmetric orientations
    for i in range(1, num_sites):
        asymmetric_orientation = "<" * i + ">" * (num_sites - i)
        orientations.append(asymmetric_orientation)

    # Add reverse of symmetric orientation
    orientations.append("<" * num_sites)
    
    # Add reverse of asymmetric orientations
    for i in range(1, num_sites):
        asymmetric_orientation_reverse = ">" * (num_sites - i) + "<" * i
        orientations.append(asymmetric_orientation_reverse)
    
    return orientations


################################################################################
# main
################################################################################

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option(
        "--input-tsv-file",
        dest="input_tsv_file",
        default="/home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv",
        help="Specify path to the file with coordinates of CTCF-binding sites in the tested genome",
    )
    parser.add_option(
        "--num-inserts",
        dest="num_inserts",
        default=6,
        type="int",
        help="Specify number of inserts (core+flanks) to be inserted into a background sequence.",
    )
    parser.add_option(
        "--flank-length",
        dest="flank_length",
        default=30,
        type="int",
        help="Specify range of right and left flank to be tested",
    )
    parser.add_option(
        "--flank-spacer-sum",
        dest="flank_spacer_sum",
        default=90,
        type="int",
        help="Specify sum of flank and spacer so that distances between CTCFs binding sites are kept constant. \
        2xflank-spacer-sum=distance between two consecutive CTCFs.",
    )
    parser.add_option(
        "--backgrounds-indices",
        dest="backgrounds_indices",
        default="0,1,2,3,4,5,6,7,8,9",
        type="string",
        help="Specify number of background sequences that CTCFs will be inserted into",
    )  
    parser.add_option(
        "--num-ctcf-sites",
        dest="num_ctcf_sites",
        default=100,
        type="int",
        help="Specify number of CTCF sites to be tested.",
    )
    parser.add_option(
        "--output-filename",
        dest="output_filename",
        default="out.tsv",
        help="Filename for output",
    )    
    
    (options, args) = parser.parse_args()
    
    background_indices_list = [int(index) for index in options.backgrounds_indices.split(",")]
    orient_list = generate_orientations(int(options.num_inserts))
    flank_start, flank_end = options.flank_length, options.flank_length

    CTCF_df = pd.read_csv(options.input_tsv_file, sep="\t")
    if "Unnamed: 0" in CTCF_df.columns:
        CTCF_df = CTCF_df.drop(columns=["Unnamed: 0"])

    CTCF_df = CTCF_df.sort_values(by="insertion_SCD", ascending=False)[:options.num_ctcf_sites]
    CTCF_df = CTCF_df.reset_index(drop=True)

    # adding orientation
    CTCF_df_with_orientation = add_orientation(
        CTCF_df,
        orientation_strings=orient_list,
        all_permutations=False,
    )

    # adding background index
    CTCF_df_with_background = add_background(
        CTCF_df_with_orientation, 
        background_indices_list
        )

    # adding flank and spacer
    CTCF_df_with_flanks_spacers = add_diff_flanks_and_const_spacer(
        CTCF_df_with_background, 
        flank_start, 
        flank_end, 
        options.flank_spacer_sum
        )

    CTCF_df_with_flanks_spacers.to_csv(
            options.output_filename, sep="\t", index=False, header=True
        )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
