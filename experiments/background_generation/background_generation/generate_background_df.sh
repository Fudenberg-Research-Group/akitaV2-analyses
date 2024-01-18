#!/bin/bash

# Conda env activation
eval "$(conda shell.bash hook)"
conda activate basenji_py3.9_tf2.15

# Define variables
genome_fasta="/project/fudenber_735/genomes/mm10/mm10.fa"
seq_bed_file="/project/fudenber_735/tensorflow_models/akita/v2/data/mm10/sequences.bed"
output_filename="./input_data/50seqs_GCuniform_maxSCD35.tsv"
SCD_threshold="35"
shuffle_parameter="8"
ctcf_detection_threshold="8"
mutation_method="permute_whole_seq"
num_backgrounds=50
mode="uniform"

# Check if required files exist
if [ ! -f "$genome_fasta" ]; then
    echo "Error: genome FASTA file does not exist."
    exit 1
fi

if [ ! -f "$seq_bed_file" ]; then
    echo "Error: seq_bed_file does not exist."
    exit 1
fi

# Run the command
python generate_background_df.py \
    -f "$genome_fasta" \
    -seq_bed_file "$seq_bed_file" \
    --output_filename "$output_filename" \
    --SCD_threshold "$SCD_threshold" \
    --shuffle_parameter "$shuffle_parameter" \
    --ctcf_detection_threshold "$ctcf_detection_threshold" \
    --mutation_method "$mutation_method" \
    --num_backgrounds "$num_backgrounds" \
    --mode "$mode" \

if [ $? -ne 0 ]; then
    echo "Error: failed to run generate_background_df.py."
    exit 1
fi

echo "All done"