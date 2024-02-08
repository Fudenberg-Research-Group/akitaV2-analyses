#!/usr/bin/env python

# Copyright 2017 Calico LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =========================================================================

###################################################

"""
generate_single_insertion_df.py

input tsv table columns:
chrom | start | end | strand

This way one row represents a single experiment.

The script requires the following input:
- orientation string
- flank range
- desired sum of the length of (flank + spacer)
- (optional) number of background sequences

"""

################################################################################
# imports
################################################################################

from optparse import OptionParser
import pandas as pd

from akita_utils.tsv_gen_utils import (
    add_orientation,
    add_background,
    add_diff_flanks_and_const_spacer,
)

import sys
sys.path.insert(0, "/home1/smaruj/akitaX1-analyses/utils/")
from analysis_utils import split_by_percentile_groups

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
    
    flank_start, flank_end = [int(flank) for flank in options.flank_range.split(",")]
    background_indices_list = [int(index) for index in options.backgrounds_indices.split(",")]
    orient_list = options.orientation_string.split(",")

    # reading the tsv with CTCF sites coordinates
    input_df = pd.read_csv(options.input_tsv_file, sep="\t")
    if "Unnamed: 0" in input_df.columns:
        input_df = input_df.drop(columns=["Unnamed: 0"])

    # spliting sites by percentiles wrt insertion_SCD
    input_df = split_by_percentile_groups(input_df, column_to_split="insertion_SCD", num_classes=5, 
                                   upper_percentile=100, lower_percentile=0, 
                                   category_colname="insSCD_group")
    # sampling 100 sites from high, medium, and low group
    high = input_df[input_df["insSCD_group"] == "Group_4"].sample(n=100)
    medium = input_df[input_df["insSCD_group"] == "Group_2"].sample(n=100)
    low = input_df[input_df["insSCD_group"] == "Group_0"].sample(n=100)
    
    # creating final tsv table
    concat_df = pd.concat([high, medium, low])
    concat_df = concat_df.reset_index(drop=True)
    
    combined_df = concat_df.merge(right=concat_df, how="cross", suffixes=("_core", "_flank"))
    
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
        options.flank_spacer_sum
        )

    # adding background index
    for background_index in background_indices_list:
        
        combined_df_with_background = add_background(
            combined_df_with_flanks_spacers, 
            [background_index]
            )
        
        combined_df_with_background.to_csv(
                options.output_filename.split(".")[0] + "bg_" + str(background_index) + "." + options.output_filename.split(".")[1], sep="\t", index=False
            )

################################################################################
# __main__
################################################################################
if __name__ == "__main__":
    main()