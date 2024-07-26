# Description:
# This script processes a TSV file of CTCF-binding site coordinates by adding orientation, background index,
# flanking regions, and varying spacer lengths. It generates a processed TSV file with additional columns
# for these features, allowing for more detailed analysis of CTCF-binding sites.
#
# Inputs:
# - --input-tsv-file: Path to the TSV file with CTCF-binding site coordinates [Default: "/home1/smaruj/akitaX1-analyses/input_data/select_top20percent/output/CTCFs_jaspar_filtered_mm10_top20percent.tsv"].
# - --orientation-string: Orientation string or list of orientation strings to be tested [Default: ">>"].
# - --flank-length: Length of the flanking regions to be tested [Default: 30].
# - --space-range: Range of spacing values to be tested [Default: "1,1000"].
# - --num_log-intervals: Number of intervals to divide the spacing range into [Default: 400].
# - --backgrounds-indices: Indices of background sequences for CTCF insertion [Default: "0,1,2,3,4,5,6,7,8,9"].
# - --output-filename: Filename for the processed output TSV file [Default: "out.tsv"].
# - --all-permutations: Boolean flag to test all possible permutations of the provided orientation string [Default: False].
#
# Example command-line usage:
# python generate_spacing_insertion_df.py --input-tsv-file input.tsv --orientation-string ">>" --flank-length 30 --space-range "1,1000" --num_log-intervals 400 --backgrounds-indices "0,1,2,3,4,5,6,7,8,9" --output-filename processed_output.tsv --all-permutations

from optparse import OptionParser
import pandas as pd
import numpy as np
from cooltools.lib import numutils

from akita_utils.tsv_utils import (
    add_orientation,
    add_background,
    add_const_flank_and_diff_spacer,
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
        default="/home1/smaruj/akitaX1-analyses/input_data/select_top20percent/output/CTCFs_jaspar_filtered_mm10_top20percent.tsv",
        help="Specify path to the file with coordinates of CTCF-binding sites in the tested genome",
    )
    parser.add_option(
        "--orientation-string",
        dest="orientation_string",
        default=">>",
        type="string",
        help="Specify orientation string - one string that will be tested for each CTCF or a list of orientation strings",
    )
    parser.add_option(
        "--flank-length",
        dest="flank_length",
        default=30,
        type="int",
        help="Specify range of right and left flank to be tested",
    )
    parser.add_option(
        "--space-range",
        dest="space_range",
        default="1,1000",
        type="string",
        help="Specify range of spacing to be tested",
    )
    parser.add_option(
        "--num_log-intervals",
        dest="log_space_range",
        default=400,
        type="int",
        help="Specify number of intervals to divide the space-range into",
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

    flank_length = options.flank_length
    background_indices_list = [
        int(index) for index in options.backgrounds_indices.split(",")
    ]
    orient_list = options.orientation_string.split(",")

    spacing_start, spacing_end = [
        int(num) for num in options.space_range.split(",")
    ]

    spacing_list = list(
        np.unique(
            numutils.logbins(
                lo=spacing_start,
                hi=spacing_end,
                N=options.log_space_range,
                version=2,
            )
            - 1
        )
    )

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
    CTCF_df_with_flanks_spacers = add_const_flank_and_diff_spacer(
        CTCF_df_with_background, flank_length, spacing_list
    )

    CTCF_df_with_flanks_spacers.to_csv(
        options.output_filename, sep="\t", index=False
    )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
