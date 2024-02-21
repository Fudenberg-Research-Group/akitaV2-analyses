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
import random
import pysam

from akita_utils.tsv_gen_utils import (
    add_orientation,
    add_background,
    add_diff_flanks_and_const_spacer,
)
from akita_utils.dna_utils import dna_rc

################################################################################
# main
################################################################################

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option(
        "-f",
        dest="genome_fasta",
        default="/project/fudenber_735/genomes/mm10/mm10.fa" ,
        help="Genome FASTA for sequences [Default: %default]",
    )
    parser.add_option(
        "--input-tsv-file",
        dest="input_tsv_file",
        default="/home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv",
        help="Specify path to the file with coordinates of CTCF-binding sites in the tested genome",
    )
    parser.add_option(
        "--lower-df-range",
        dest="lower_df_range",
        default=1500,
        type="int",
        help="Specify part of df you're interested in",
    )
    parser.add_option(
        "--upper-df-range",
        dest="upper_df_range",
        default=0,
        type="int",
        help="Specify part of df you're interested in",
    )
    parser.add_option(
        "--orientation-string",
        dest="orientation_string",
        default=">",
        type="string",
        help="Specify orientation string - one string that will be tested for each CTCF or a list of orientation strings",
    )
    parser.add_option(
        "--flank-length",
        dest="flank_length",
        default=15,
        type="int",
        help="Specify length of flanking sequences",
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

    background_indices_list = [int(index) for index in options.backgrounds_indices.split(",")]
    orient_list = options.orientation_string.split(",")
    flank_length = options.flank_length
    
    CTCF_df = pd.read_csv(options.input_tsv_file, sep="\t")
    if "Unnamed: 0" in CTCF_df.columns:
        CTCF_df = CTCF_df.drop(columns=["Unnamed: 0"])
    CTCF_df = CTCF_df.sort_values(by="insertion_SCD", ascending=False)[options.upper_df_range:options.lower_df_range]

    # open genome FASTA
    genome_open = pysam.Fastafile(
        options.genome_fasta
    )
    
    # Initialize a list to store experiment data
    experiment_data = []
    
    # Iterate over each site
    for index, row in CTCF_df.iterrows():
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]
        insertion_SCD = row["insertion_SCD"]
        disruption_SCD = row["disruption_SCD"]
        
        # Adjust start and end positions to include flanking regions
        start_adj = start - flank_length
        end_adj = end + flank_length
    
        # Fetch the sequence from the genome
        sequence = genome_open.fetch(chrom, start_adj, end_adj).upper()
    
        # If the site is on the negative strand, reverse complement the sequence
        if strand == '-':
            sequence = dna_rc(sequence).upper()

        # Iterate over each position in the sequence
        for i in range(len(sequence)):
            original_nucleotide = sequence[i]
            # Mutate the nucleotide to any of the other three nucleotides
            for mutated_nucleotide in ['A', 'T', 'G', 'C']:
                # Append to the experiment data list, including inherited information
                experiment_data.append([
                    chrom, start, end, strand, insertion_SCD, disruption_SCD,
                    (i - flank_length), original_nucleotide, mutated_nucleotide
                ])
        
    # Convert the experiment data list to a DataFrame
    columns = ["chrom", "start", "end", "strand", "insertion_SCD", "disruption_SCD",
               "position", "original_nucleotide", "mutated_nucleotide"]
    experiment_df = pd.DataFrame(experiment_data, columns=columns)

    # Adding insertion-specific columns
    experiment_df_with_orientation = add_orientation(
        experiment_df,
        orientation_strings=orient_list,
        all_permutations=False,
    )

    # adding background index
    experiment_df_with_background = add_background(experiment_df_with_orientation, 
                                                background_indices_list
                                                )

    # adding flank and spacer
    experiment_df_with_flanks_spacers = add_diff_flanks_and_const_spacer(
        experiment_df_with_background, 
        flank_length, 
        flank_length, 
        options.flank_spacer_sum
        )

    # Write to a new TSV file
    experiment_df_with_flanks_spacers.to_csv(options.output_filename, sep="\t", index=False)


################################################################################
# __main__
################################################################################
if __name__ == "__main__":
    main()
