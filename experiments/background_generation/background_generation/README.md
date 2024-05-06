# Background Generation

This directory facilitates the generation of background genomic data for analysis.

## Directories:

### 1. input_data/
Contains a TSV table with 50 genomic windows with uniformly distributed GC content.

## Files:

- **generate_background_df.py**: Given parameters (like a bed file of sequences, shuffle parameter, signal threshold), it prepares a dataframe with genomic windows, calculates their GC content, and appends the dataframe with other specified parameters.
- **generate_background_df.sh**: Bash file running the above script under the slurm system (sbatch).
- **generate_background_sequences.py**: Given a TSV table specifying genomic windows and mutation parameters, it creates flat background sequences.
- **multiGPU-generate_background_sequences.sh**: Python script automating the above script on multiple GPUs.
- **generate_background_sequences.sh**: Script to run the above multi-GPU python script under the slurm system.