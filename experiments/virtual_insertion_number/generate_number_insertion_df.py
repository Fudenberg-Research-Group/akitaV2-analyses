from optparse import OptionParser
import pandas as pd

from akita_utils.tsv_utils import (
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
        "--max-number-sites",
        dest="max_number_sites",
        default=20,
        type="int",
        help="Specify maximum number of sites inserted",
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
        "--output-filename",
        dest="output_filename",
        default="out.tsv",
        help="Filename for output",
    )    

    (options, args) = parser.parse_args()
    
    background_indices_list = [int(index) for index in options.backgrounds_indices.split(",")]
    orient_list = [x*options.orientation_string for x in range(1, options.max_number_sites+1)]
    
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
        options.flank_length, 
        options.flank_length, 
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

