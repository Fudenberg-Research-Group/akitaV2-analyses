# Description:
# This script processes a TSV file containing coordinates of CTCF-binding sites in a specified genome.
# It adds background sequences, orientations, and spacers to generate two scenarios: boundary and dot.
# The processed data is saved to separate output files for each scenario.
#
# Inputs:
# - input_tsv_file: Path to the input TSV file with coordinates of CTCF-binding sites.
# - flank_length: Length of flanking sequences to be added around CTCF sites.
# - boundary_orientation_string: Orientation string for the boundary scenario.
# - boundary_spacer: Spacer length for the boundary scenario.
# - dot_orientation_string: Orientation string for the dot scenario.
# - dot_spacer: Spacer length for the dot scenario.
# - backgrounds_indices: Comma-separated list of background sequence indices for CTCF insertion.
# - boundary_output_filename: Filename for the boundary scenario output.
# - dot_output_filename: Filename for the dot scenario output.
#
# Outputs:
# - Two TSV files containing processed CTCF-binding site data with added background sequences, orientations, and spacers:
#   1. boundary_output_filename: Output file for the boundary scenario.
#   2. dot_output_filename: Output file for the dot scenario.
#
# Example command-line usage:
# python generate_tsv_dot_boundary_scenario.py --input-tsv-file CTCFs.tsv --flank-length 30 --boundary-orientation-string "<>"
# --boundary-spacer 60 --dot-orientation-string "><" --dot-spacer 199970 --backgrounds-indices "0,1,2,3,4,5,6,7,8,9"
# --boundary-output-filename boundary_out.tsv --dot-output-filename dot_out.tsv

from optparse import OptionParser
import pandas as pd

from akita_utils.df_utils import (
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
        default="/home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv",
        help="Specify path to the file with coordinates of CTCF-binding sites in the tested genome",
    )
    parser.add_option(
        "--flank-lenght",
        dest="flank_lenght",
        default=30,
        type="int",
        help="Specify lenght of flanking sequences",
    )
    parser.add_option(
        "--boundary-orientation-string",
        dest="boundary_orientation_string",
        default="<>",
        type="string",
        help="Specify orientation string for the boundary scenario",
    )
    parser.add_option(
        "--boundary-spacer",
        dest="boundary_spacer",
        default=60,
        type="int",
        help="Specify spacer for the boundary scenario",
    )
    parser.add_option(
        "--dot-orientation-string",
        dest="dot_orientation_string",
        default="><",
        type="string",
        help="Specify orientation string for the dot scenario",
    )
    parser.add_option(
        "--dot-spacer",
        dest="dot_spacer",
        default=199970,
        type="int",
        help="Specify spacer for the dot scenario",
    )
    parser.add_option(
        "--backgrounds-indices",
        dest="backgrounds_indices",
        default="0,1,2,3,4,5,6,7,8,9",
        type="string",
        help="Specify number of background sequences that CTCFs will be inserted into",
    )
    parser.add_option(
        "--boundary-output-filename",
        dest="boundary_output_filename",
        default="boundary_out.tsv",
        help="Filename for boundary output",
    )
    parser.add_option(
        "--dot-output-filename",
        dest="dot_output_filename",
        default="dot_out.tsv",
        help="Filename for dot output",
    )

    (options, args) = parser.parse_args()

    flank_length = options.flank_lenght
    background_indices_list = [
        int(index) for index in options.backgrounds_indices.split(",")
    ]

    boundary_orient_list = [options.boundary_orientation_string]
    dot_orient_list = [options.dot_orientation_string]

    boundary_spacing_list = [options.boundary_spacer]
    dot_spacing_list = [options.dot_spacer]

    CTCF_df = pd.read_csv(options.input_tsv_file, sep="\t")
    nr_sites = len(CTCF_df)

    # adding background index
    CTCF_df_with_background = add_background(CTCF_df, background_indices_list)

    exp_id = [i for i in range(nr_sites * len(background_indices_list))]
    CTCF_df_with_background["exp_id"] = exp_id

    # adding orientation for boundary
    boundary_CTCF_df_with_orientation = add_orientation(
        CTCF_df_with_background.copy(),
        orientation_strings=boundary_orient_list,
        all_permutations=False,
    )

    # adding orientation for dots
    dot_CTCF_df_with_orientation = add_orientation(
        CTCF_df_with_background.copy(),
        orientation_strings=dot_orient_list,
        all_permutations=False,
    )

    # adding flank and spacer for boundary
    boundary_CTCF_df_with_flanks_spacers = add_const_flank_and_diff_spacer(
        boundary_CTCF_df_with_orientation, flank_length, boundary_spacing_list
    )

    # adding flank and spacer for dot
    dot_CTCF_df_with_flanks_spacers = add_const_flank_and_diff_spacer(
        dot_CTCF_df_with_orientation, flank_length, dot_spacing_list
    )

    # saving
    boundary_CTCF_df_with_flanks_spacers.to_csv(
        options.boundary_output_filename, sep="\t", index=False
    )

    dot_CTCF_df_with_flanks_spacers.to_csv(
        options.dot_output_filename, sep="\t", index=False
    )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
