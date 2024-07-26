# Description:
# This script processes CTCF-binding sites from a given genome and generates a TSV file with various experimental data.
# It fetches sequences from the genome around the CTCF-binding sites, mutates nucleotides, and adds experimental context such as orientation, background indices, and flank/spacer information.
#
# Inputs:
# - --input-tsv-file: Path to the TSV file with CTCF-binding sites.
# - -f, --genome-fasta: Path to the genome FASTA file [Default: "/project/fudenber_735/genomes/mm10/mm10.fa"].
# - --lower-df-range: Lower range for filtering CTCF sites by SCD [Default: 1500].
# - --upper-df-range: Upper range for filtering CTCF sites by SCD [Default: 0].
# - --orientation-string: String(s) specifying orientation [Default: ">"].
# - --flank-length: Length of flanking sequences [Default: 15].
# - --flank-spacer-sum: Sum of flank and spacer lengths [Default: 30].
# - --backgrounds-indices: Indices of background sequences [Default: "0,1,2,3,4,5,6,7,8,9"].
# - --output-filename: Filename for the output TSV file [Default: "out.tsv"].
#
# Outputs:
# - A TSV file with experimental data, including orientation, background indices, and flank/spacer information.
#
# Example command-line usage:
# python generate_pairwise_mutagenesis_df.py --input-tsv-file CTCFs.tsv -f genome.fa --lower-df-range 1500 --upper-df-range 0 --orientation-string ">" --flank-length 15 --flank-spacer-sum 30 --backgrounds-indices "0,1,2,3,4" --output-filename experiment_data.tsv

from optparse import OptionParser
import pandas as pd
import random
import pysam

from akita_utils.tsv_utils import (
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
        default="/project/fudenber_735/genomes/mm10/mm10.fa",
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

    background_indices_list = [
        int(index) for index in options.backgrounds_indices.split(",")
    ]
    orient_list = options.orientation_string.split(",")
    flank_length = options.flank_length

    CTCF_df = pd.read_csv(options.input_tsv_file, sep="\t")
    if "Unnamed: 0" in CTCF_df.columns:
        CTCF_df = CTCF_df.drop(columns=["Unnamed: 0"])
    CTCF_df = CTCF_df.sort_values(by="insertion_SCD", ascending=False)[
        options.upper_df_range : options.lower_df_range
    ]

    # open genome FASTA
    genome_open = pysam.Fastafile(options.genome_fasta)

    # Function to generate a random nucleotide different from the original
    def mutate_nucleotide(original):
        nucleotides = ["A", "T", "G", "C"]
        nucleotides.remove(original)
        return random.choice(nucleotides)

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
        if strand == "-":
            sequence = dna_rc(sequence).upper()

        for i in range(len(sequence)):
            for j in range(
                i + 1, len(sequence)
            ):  # Pairwise combinations without repetition
                pos1 = i - flank_length
                pos2 = j - flank_length
                original_nuc1 = sequence[i]
                original_nuc2 = sequence[j]
                mutated_nuc1 = mutate_nucleotide(original_nuc1)
                mutated_nuc2 = mutate_nucleotide(original_nuc2)
                # Append to the experiment data list, including inherited information
                experiment_data.append(
                    [
                        chrom,
                        start,
                        end,
                        strand,
                        insertion_SCD,
                        disruption_SCD,
                        pos1,
                        pos2,
                        original_nuc1,
                        mutated_nuc1,
                        original_nuc2,
                        mutated_nuc2,
                    ]
                )

    # Convert the experiment data list to a DataFrame
    columns = [
        "chrom",
        "start",
        "end",
        "strand",
        "insertion_SCD",
        "disruption_SCD",
        "pos1",
        "pos2",
        "OriginalNuc1",
        "MutatedNuc1",
        "OriginalNuc2",
        "MutatedNuc2",
    ]
    experiment_df = pd.DataFrame(experiment_data, columns=columns)

    # Adding insertion-specific columns
    experiment_df_with_orientation = add_orientation(
        experiment_df,
        orientation_strings=orient_list,
        all_permutations=False,
    )

    # adding background index
    experiment_df_with_background = add_background(
        experiment_df_with_orientation, background_indices_list
    )

    # adding flank and spacer
    experiment_df_with_flanks_spacers = add_diff_flanks_and_const_spacer(
        experiment_df_with_background,
        flank_length,
        flank_length,
        options.flank_spacer_sum,
    )

    # Write to a new TSV file
    experiment_df_with_flanks_spacers.to_csv(
        options.output_filename, sep="\t", index=False
    )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
