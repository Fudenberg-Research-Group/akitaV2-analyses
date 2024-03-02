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
generate_insulation_offset_assymetry_df.py

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
    add_diff_flanks_and_const_spacer,
)


################################################################################
# helper function
################################################################################


def generate_orientations_for_multiple_sites(site_numbers):
    """
    Generate orientations for different numbers of CTCF sites.

    Parameters
    ------------
    site_numbers : list
        List of numbers of CTCF sites for which orientations need to be generated.

    Returns
    ---------
    orientations : list
        List containing orientations for all specified numbers of sites.
    """
    orientations = []

    for num_sites in site_numbers:
        
        # Generate symmetric orientations
        symmetric_orientation_1 = "<" * (num_sites // 2) + ">" * (num_sites // 2)
        symmetric_orientation_2 = ">" * (num_sites // 2) + "<" * (num_sites // 2)
        orientations.extend([symmetric_orientation_1, symmetric_orientation_2])
        
        # Generate asymmetric orientations
        asymmetric_orientation_1 = "<" * num_sites
        asymmetric_orientation_2 = ">" * num_sites
        orientations.extend([asymmetric_orientation_1, asymmetric_orientation_2])

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
        default="2,4,6,8,10,12",
        type="str",
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
    site_numbers = [int(index) for index in options.num_inserts.split(",")]
    orient_list = generate_orientations_for_multiple_sites(site_numbers)
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

