# Description:
# This script processes a TSV file containing CTCF binding site coordinates by adding orientation strings, 
# background indices, and flank/spacer information. The processed data is saved to a new TSV file.
#
# Inputs:
# - --input-tsv-file: Path to the input TSV file with CTCF binding site coordinates [Default: "/home1/smaruj/akitaX1-analyses/input_data/preprocess_CTCFs/output/CTCFs_jaspar_filtered_mm10.tsv"].
# - --orientation-string: String or comma-separated list of strings specifying orientation for CTCF sites [Default: ">"].
# - --flank-range: Range of flanks to be tested (left and right) [Default: "30"].
# - --flank-spacer-sum: Sum of flank and spacer distances to keep spacing constant [Default: 30].
# - --backgrounds-indices: Comma-separated list of indices for background sequences [Default: "0,1,2,3,4,5,6,7,8,9"].
# - --output-filename: Filename for the output TSV file [Default: "out.tsv"].
#
# Outputs:
# - Processed TSV file with additional columns for orientation, background index, and flank/spacer information.
#
# Example command-line usage:
# python generate_single_insertion_df.py --input-tsv-file input_data.tsv --orientation-string ">" --flank-range 30 --flank-spacer-sum 30 --backgrounds-indices "0,1,2" --output-filename processed_ctcf_data.tsv

from optparse import OptionParser
import pandas as pd

from akita_utils.df_utils import (
    add_orientation,
    add_background,
    add_diff_flanks_and_const_spacer,
)


################################################################################
# main
################################################################################

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option(
        "--input-tsv-file",
        dest="input_tsv_file",
        default="/home1/smaruj/akitaX1-analyses/input_data/preprocess_CTCFs/output/CTCFs_jaspar_filtered_mm10.tsv",
        help="Specify path to the file with coordinates of CTCF-binding sites in the tested genome",
    )
    parser.add_option(
        "--orientation-string",
        dest="orientation_string",
        default=">",
        type="string",
        help="Specify orientation string - one string that will be tested for each CTCF or a list of orientation strings",
    )
    parser.add_option(
        "--flank-range",
        dest="flank_range",
        default="30",
        type="string",
        help="Specify range of right and left flank to be tested",
    )
    parser.add_option(
        "--flank-spacer-sum",
        dest="flank_spacer_sum",
        default=30,
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
        "--output-filename",
        dest="output_filename",
        default="out.tsv",
        help="Filename for output",
    )    

    (options, args) = parser.parse_args()
    
    flank_start, flank_end = [int(options.flank_range), int(options.flank_range)]
    background_indices_list = [int(index) for index in options.backgrounds_indices.split(",")]
    orient_list = options.orientation_string.split(",")
    
    CTCF_df = pd.read_csv(options.input_tsv_file, sep="\t")

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
            options.output_filename, sep="\t", index=False
        )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
