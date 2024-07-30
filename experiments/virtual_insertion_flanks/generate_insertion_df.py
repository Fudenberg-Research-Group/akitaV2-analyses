# Description:
# This script processes a TSV file containing CTCF-binding site coordinates by adding orientation, background indices, and flank/spacer information.
# The resulting data is saved to a new TSV file.
#
# Inputs:
# - --input-tsv-file: Path to the input TSV file with CTCF-binding site coordinates [Default: /home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv].
# - --orientation-string: Orientation string(s) to be tested for each CTCF [Default: ">"].
# - --flank-range: Range of right and left flank to be tested, specified as "left_flank,right_flank" [Default: "0,35"].
# - --flank-spacer-sum: Sum of flank and spacer distances to keep distances between CTCFs constant [Default: 90].
# - --backgrounds-indices: Comma-separated indices of background sequences [Default: "0,1,2,3,4,5,6,7,8,9"].
# - --output-filename: Filename for the output TSV file [Default: "out.tsv"].
# - --all-permutations: Boolean flag to test all possible permutations of the provided orientation string(s) [Default: False].
#
# Outputs:
# - TSV file with additional columns for orientation, background index, and flank/spacer information.
#
# Example command-line usage:
# python generate_insertion_df.py --input-tsv-file input.tsv --orientation-string "><" --flank-range "0,20" --flank-spacer-sum 100 --backgrounds-indices "0,1,2" --output-filename processed_ctcfs.tsv --all-permutations


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
        default="/home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv",
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
        default="0,35",
        type="string",
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
        "--output-filename",
        dest="output_filename",
        default="out.tsv",
        help="Filename for output",
    )
    parser.add_option(
        "--all-permutations",
        dest="all_permutations",
        default=False,
        action="store_true",
        help="Test all possible permutations of N = length of provided orientation_string",
    )

    (options, args) = parser.parse_args()

    flank_start, flank_end = [
        int(flank) for flank in options.flank_range.split(",")
    ]
    background_indices_list = [
        int(index) for index in options.backgrounds_indices.split(",")
    ]
    orient_list = options.orientation_string.split(",")

    CTCF_df = pd.read_csv(options.input_tsv_file, sep="\t")

    # adding orientation
    CTCF_df_with_orientation = add_orientation(
        CTCF_df,
        orientation_strings=orient_list,
        all_permutations=options.all_permutations,
    )

    # adding background index
    CTCF_df_with_background = add_background(
        CTCF_df_with_orientation, background_indices_list
    )

    # adding flank and spacer
    CTCF_df_with_flanks_spacers = add_diff_flanks_and_const_spacer(
        CTCF_df_with_background,
        flank_start,
        flank_end,
        options.flank_spacer_sum,
    )

    CTCF_df_with_flanks_spacers.to_csv(
        options.output_filename, sep="\t", index=False
    )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
