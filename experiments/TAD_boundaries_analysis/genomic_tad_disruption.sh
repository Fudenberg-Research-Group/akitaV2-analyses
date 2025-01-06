#!/bin/bash

# Conda env activation
eval "$(conda shell.bash hook)"
conda activate basenji_py3.9_tf2.15

# Parse command line arguments, replace with custom parameters
genome_fasta="/project/fudenber_735/genomes/mm10/mm10.fa" 
models_dir="/project/fudenber_735/tensorflow_models/akita/v2/models"
tsv_file="/home1/smaruj/akitaV2-analyses/experiments/TAD_boundaries_analysis/input_data/boundaries_disruption_1.tsv"
out_dir="/scratch1/smaruj/tad_disruption/model0_1" 
models="0" # this is a string with space seperated model numbers examples "0 1 2" or "4 5 6" 
batch_size=8 
max_proc=30
processes=30
stats="SCD"
time="0-02:30:00" 
# constraint="[xeon-6130|xeon-2640v4]"

# Check if genome_fasta file exists
if [ ! -f "$genome_fasta" ]; then
    echo "Genome fasta file does not exist."
    exit 1
fi

# Check if models_dir directory exists
if [ ! -d "$models_dir" ]; then
    echo "Models directory does not exist."
    exit 1
fi

# Check if tsv_file file exists
if [ ! -f "$tsv_file" ]; then
    echo "TSV file does not exist."
    exit 1
fi

for model in $models

do 
    # infering params_file and model_file paths
    this_models_dir="${models_dir}";
    this_models_dir+="/f";
    this_models_dir+="${model}";
    this_models_dir+="c0/train/";
    
    params_file="${this_models_dir}";
    params_file+="params.json";
    
    model_file="${this_models_dir}";
    model_file+="model1_best.h5"; #mouse
    
    # changing outdir for this model
    this_out_dir="${out_dir}";
    this_out_dir+="_m";
    this_out_dir+="${model}";

    # changing name of the jobs to be submitted
    name="D_m";
    name+="${model}";
    
    # running multiGPU script
    python multiGPU-genomic_tad_disruption.py "${params_file}" "${model_file}" "${tsv_file}" -f "${genome_fasta}" -o "${this_out_dir}" --stats "${stats}" --batch-size "${batch_size}" -p "${processes}" --max_proc "${max_proc}" --name "${name}" --time "${time}"
    sleep 15
done

echo All jobs submitted

