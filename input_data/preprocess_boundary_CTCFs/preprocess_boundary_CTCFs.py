"""
This script processes genomic data to identify CTCF binding sites that overlap
with genomic boundaries and are filtered based on various criteria. The script
uses the Akita deep learning model to predict genome folding from DNA sequence.
It filters CTCF motifs, repetitive elements, and genomic boundaries, then 
identifies overlaps between CTCF binding sites and boundaries. The final output
is a TSV file with the filtered and overlapping CTCF binding sites.

Usage:
    python script.py [options] <params_file> <vcf_file>

Options:
    --model-params-file         Path to the model parameters JSON file.
                                [Default: /project/fudenber_735/tensorflow_models/akita/v2/models/f0c0/train/params.json]
    --jaspar-file               Path to the JASPAR file with CTCF sites coordinates.
                                [Default: /project/fudenber_735/motifs/mm10/jaspar/MA0139.1.tsv.gz]
    --ctcf-filter-expand-window Window size for CTCF filtering. [Default: 60]
    --rmsk-file                 Path to the RMSK file. [Default: /project/fudenber_735/genomes/mm10/database/rmsk.txt.gz]
    --rmsk-filter-expand-window Window size for RMSK filtering. [Default: 20]
    --chrom-sizes-file          Path to the chromosome sizes file. [Default: /project/fudenber_735/genomes/mm10/mm10.chrom.sizes.reduced]
    --boundary-file             Path to the boundary file. [Default: /project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/features/bonev2017.HiC_ES.mm10.mapq_30.1000.window_200000.insulation]
    --boundary-strength-thresh  Threshold for boundary strength. [Default: 0.25]
    --boundary-insulation-thresh Threshold for boundary insulation score. [Default: 0.00]
    --output-tsv-path           Output path for the filtered CTCF sites TSV file. [Default: ./output/CTCFs_jaspar_filtered_mm10.tsv]
    --autosomes-only            Drop the sex chromosomes and mitochondrial DNA. [Default: True]

The script performs the following steps:
1. Parses command-line options and arguments.
2. Loads the model parameters to obtain the sequence length.
3. Loads JASPAR CTCF motifs and filters them based on chromosome ID if specified.
4. Reads the RMSK file for repetitive elements.
5. Loads and filters boundaries based on strength and insulation thresholds.
6. Identifies overlaps between filtered boundaries and CTCF motifs.
7. Filters overlapping CTCF sites by the number of overlaps with other CTCF sites and RMSK elements.
8. Adds a sequence ID to the filtered CTCF sites.
9. Saves the final filtered CTCF sites to a TSV file.
"""


from optparse import OptionParser
import json
import bioframe as bf
import pandas as pd
import os
from akita_utils.df_utils import (
    filter_by_chrmlen,
    filter_by_overlap_num,
    filter_by_chromID,
)
from akita_utils.format_io import read_rmsk


def main():
    usage = "usage: %prog [options] <params_file> <vcf_file>"
    parser = OptionParser(usage)

    parser.add_option(
        "--model-params-file",
        dest="model_params_file",
        default="/project/fudenber_735/tensorflow_models/akita/v2/models/f0c0/train/params.json",
        help="Parameters of model to be used[Default: %default]",
    )
    parser.add_option(
        "--jaspar-file",
        dest="jaspar_file",
        default="/project/fudenber_735/motifs/mm10/jaspar/MA0139.1.tsv.gz",
        help="Jaspar file with ctcf sites coordinates [Default: %default]",
    )
    parser.add_option(
        "--ctcf-filter-expand-window",
        dest="ctcf_filter_expand_window",
        default=60,
        type=int,
        help="window size for the ctcf-filtering [Default: %default]",
    )
    parser.add_option(
        "--rmsk-file",
        dest="rmsk_file",
        default="/project/fudenber_735/genomes/mm10/database/rmsk.txt.gz",
        help=" [Default: %default]",
    )
    parser.add_option(
        "--rmsk-filter-expand-window",
        dest="rmsk_filter_expand_window",
        default=20,
        type=int,
        help="window size for the rmsk-filtering [Default: %default]",
    )
    parser.add_option(
        "--chrom-sizes-file",
        dest="chrom_sizes_file",
        default="/project/fudenber_735/genomes/mm10/mm10.chrom.sizes.reduced",
        help=" [Default: %default]",
    )
    parser.add_option(
        "--boundary-file",
        dest="boundary_file",
        default="/project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/features/bonev2017.HiC_ES.mm10.mapq_30.1000.window_200000.insulation",
        help=" [Default: %default]",
    )
    parser.add_option(
        "--boundary-strength-thresh",
        dest="boundary_strength_thresh",
        default=0.25,
        type=float,
        help="threshold on boundary strengths [Default: %default]",
    )
    parser.add_option(
        "--boundary-insulation-thresh",
        dest="boundary_insulation_thresh",
        default=0.00,
        type=float,
        help="threshold on boundary insulation score [Default: %default]",
    )
    parser.add_option(
        "--output-tsv-path",
        dest="output_tsv_path",
        default="./output/CTCFs_jaspar_filtered_mm10.tsv",
        type="str",
        help="Output path [Default: %default]",
    )
    parser.add_option(
        "--autosomes-only",
        dest="autosomes_only",
        default=True,
        action="store_true",
        help="Drop the sex chromosomes and mitochondrial DNA [Default: %default]",
    )

    (options, args) = parser.parse_args()

    if options.autosomes_only:
        chromID_to_drop = ["chrX", "chrY", "chrM"]

    if os.path.exists(options.output_tsv_path) is True:
        raise ValueError("boundary file already exists!")

    # get model seq_length
    with open(options.model_params_file) as params_open:
        params_model = json.load(params_open)["model"]
        seq_length = params_model["seq_length"]
    if seq_length != 1310720:
        raise Warning("potential incompatibilities with AkitaV2 seq_length")

    # load jaspar CTCF motifs
    jaspar_df = bf.read_table(options.jaspar_file, schema="jaspar", skiprows=1)
    if options.autosomes_only:
        jaspar_df = filter_by_chromID(jaspar_df, chrID_to_drop=chromID_to_drop)
    jaspar_df.reset_index(drop=True, inplace=True)

    # read rmsk file
    rmsk_df = read_rmsk(options.rmsk_file)

    # load boundaries and use standard filters for their strength
    boundaries = pd.read_csv(options.boundary_file, sep="\t")

    window_size = options.boundary_file.split("window_")[1].split(".")[0]
    boundary_key, insulation_key = (
        f"boundary_strength_{window_size}",
        f"log2_insulation_score_{window_size}",
    )

    boundaries = boundaries.iloc[
        (boundaries[boundary_key].values > options.boundary_strength_thresh)
        * (
            boundaries[insulation_key].values
            < options.boundary_insulation_thresh
        )
    ]

    if options.autosomes_only:
        boundaries = filter_by_chromID(
            boundaries, chrID_to_drop=chromID_to_drop
        )

    boundaries = filter_by_chrmlen(
        boundaries,
        options.chrom_sizes_file,
        seq_length,
    )

    boundaries.reset_index(drop=True, inplace=True)

    # overlapping CTCF df with boundaries df
    df_overlap = bf.overlap(
        boundaries, jaspar_df, suffixes=("", "_2"), return_index=False
    )

    # removing rows with no start and end info
    df_overlap = df_overlap[pd.notnull(df_overlap["start_2"])]
    df_overlap = df_overlap[pd.notnull(df_overlap["end_2"])]

    df_overlap["span"] = (
        df_overlap["start"].astype(str) + "-" + df_overlap["end"].astype(str)
    )

    df_keys = [
        "chrom",
        "start_2",
        "end_2",
        "span",
        "score_2",
        "strand_2",
        insulation_key,
        boundary_key,
    ]

    df_overlap = df_overlap[df_keys]

    # renaming
    df_overlap = df_overlap.rename(
        columns={
            "span": "boundary_span",
            "score_2": "jaspar_score",
            "start_2": "start",
            "end_2": "end",
            "strand_2": "strand",
        }
    )

    # filtering by CTCF
    filtered_df = filter_by_overlap_num(
        df_overlap,
        filter_df=jaspar_df,
        expand_window=options.ctcf_filter_expand_window,
        max_overlap_num=1,
    )

    # filtering by rmsk
    filtered_df = filter_by_overlap_num(
        filtered_df,
        rmsk_df,
        expand_window=options.rmsk_filter_expand_window,
        working_df_cols=["chrom", "start", "end"],
    )

    # adding seq_id
    filtered_df["seq_id"] = [
        seq_index for seq_index in range(len(filtered_df))
    ]

    # saving
    filtered_df.to_csv(options.output_tsv_path, sep="\t", index=False)


################################################################################
# __main__
################################################################################
if __name__ == "__main__":
    main()
