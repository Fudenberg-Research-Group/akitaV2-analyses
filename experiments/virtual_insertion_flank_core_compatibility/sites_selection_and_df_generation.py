# Description:
# This script processes a TSV file containing CTCF-binding site coordinates and generates new TSV files with various permutations of the input data. It adds orientation, flanks, spacers, and background indices to the data and saves the results into separate output files for each background index.
#
# Inputs:
# - --input-tsv-file: Path to the TSV file with coordinates of CTCF-binding sites in the genome [Default: /home1/smaruj/akitaX1-analyses/input_data/downsample_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv].
# - --orientation-string: String or list of strings specifying the orientation to be tested for each CTCF [Default: >].
# - --flank-range: Range of right and left flanks to be tested, given as "left_flank,right_flank" [Default: 30,30].
# - --flank-spacer-sum: Sum of flank and spacer distances to keep the distances between consecutive CTCFs constant [Default: 90].
# - --backgrounds-indices: Comma-separated list of background sequence indices into which CTCFs will be inserted [Default: 0,1,2,3,4,5,6,7,8,9].
# - --output-filename: Filename for the output TSV files [Default: out.tsv].
# - --all-permutations: Whether to test all possible permutations of the provided orientation string [Default: False].
#
# Outputs:
# - Multiple TSV files, each named with a background index, containing processed CTCF-binding site data.
#
# Example command-line usage:
# python sites_selection_and_df_generation.py --input-tsv-file /path/to/input.tsv --orientation-string "><" --flank-range "20,40" --flank-spacer-sum 80 --backgrounds-indices "0,1,2" --output-filename processed_data.tsv --all-permutations

from optparse import OptionParser
import pandas as pd

from akita_utils.df_utils import (
    add_orientation,
    add_background,
    add_diff_flanks_and_const_spacer,
)
from akita_utils.analysis_utils import split_by_percentile_groups


################################################################################
# main
################################################################################


def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option(
        "--input-tsv-file",
        dest="input_tsv_file",
        default="/home1/smaruj/akitaX1-analyses/input_data/downsample_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv",
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
        default="30,30",
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

    # reading the tsv with CTCF sites coordinates
    input_df = pd.read_csv(options.input_tsv_file, sep="\t")
    if "Unnamed: 0" in input_df.columns:
        input_df = input_df.drop(columns=["Unnamed: 0"])

    # spliting sites by percentiles wrt insertion_SCD
    input_df = split_by_percentile_groups(
        input_df,
        column_to_split="insertion_SCD",
        num_classes=5,
        upper_percentile=100,
        lower_percentile=0,
        category_colname="insSCD_group",
    )
    # sampling 100 sites from high, medium, and low group
    high = input_df[input_df["insSCD_group"] == "Group_4"].sample(n=100)
    medium = input_df[input_df["insSCD_group"] == "Group_2"].sample(n=100)
    low = input_df[input_df["insSCD_group"] == "Group_0"].sample(n=100)

    # creating final tsv table
    concat_df = pd.concat([high, medium, low])
    concat_df = concat_df.reset_index(drop=True)

    combined_df = concat_df.merge(
        right=concat_df, how="cross", suffixes=("_core", "_flank")
    )

    # adding orientation
    combined_df_with_orientation = add_orientation(
        combined_df,
        orientation_strings=orient_list,
        all_permutations=options.all_permutations,
    )

    # adding flank and spacer
    combined_df_with_flanks_spacers = add_diff_flanks_and_const_spacer(
        combined_df_with_orientation,
        flank_start,
        flank_end,
        options.flank_spacer_sum,
    )

    # adding background index
    for background_index in background_indices_list:
        combined_df_with_background = add_background(
            combined_df_with_flanks_spacers, [background_index]
        )

        combined_df_with_background.to_csv(
            options.output_filename.split(".")[0]
            + "bg_"
            + str(background_index)
            + "."
            + options.output_filename.split(".")[1],
            sep="\t",
            index=False,
        )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
