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
generate_spacing_insertion_df.py

input tsv table columns:
chrom | start | end | strand

This way one row represents a single experiment.

The script requires the following input:
- orientation string
- flank length
- desired sum of the length of (flank + spacer)
- (optional) number of background sequences

"""

################################################################################
# imports
################################################################################

from optparse import OptionParser
import pandas as pd
import numpy as np

from cooltools.lib import numutils

from akita_utils.tsv_gen_utils import (
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
    background_indices_list = [int(index) for index in options.backgrounds_indices.split(",")]
    orient_list = options.orientation_string.split(",")

    spacing_start, spacing_end = [int(num) for num in options.space_range.split(",")]

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
        CTCF_df_with_orientation, 
        background_indices_list
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

