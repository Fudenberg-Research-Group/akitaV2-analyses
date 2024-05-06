# Sequence Shuffling

This directory facilitates the analysis and generation of shuffled genomic sequences.

## Directories:

### 1. analysis/
Contains notebooks and plots related to sequence shuffling analysis.
- **plots/**: Directory with saved plots.
- **analysis_shuffled_vs_unshuffled.ipynb**: Notebook analyzing the insertion of a strong CTCF site into genomic sequences and the same sequences shuffled once (with shuffling_parameter=8).
- **analysis_shuffling_parameter.ipynb**: Notebook analyzing which shuffle parameter gives the flattest maps.

### 2. input_data/
Contains input data for sequence shuffling analysis.
- **shuffled_600seqs.tsv**: Hundreds (n=590) of genomic windows with shuffling parameters specified as 1, 2, 4, 8, 16, 32.
- **shuffled_600seqs_k8.tsv**: Hundreds (n=590) of genomic windows with shuffling parameter specified as 8.

## Notebooks:

- **generate_shuffled_seqs_df.py**: Notebook generating TSV with genomic windows given a set of parameters.

### For the shuffling parameter analysis:
- **generate_scores_for_shuffled_seqs.py**: Given TSV with genomic windows, this script calculates their predicted signal strength (SCD).
- **multiGPU-generate_scores_for_shuffled_seqs.py**: Runs the above script on multiple GPUs.
- **generate_scores_for_shuffled_seqs.sh**: Automates generating signal strength scores under the slurm system.

### For the shuffled vs. unshuffled comparison:
- **generate_scores_for_shuffled_insertions.py & generate_scores_for_unshuffled_insertions.py**: Given TSV with genomic windows, these scripts generate predicted signal strength scores (SCD) for insertions into respectively unshuffled and shuffled (once with k=8) genomic windows.
- **multiGPU-generate_scores_for_shuffled_insertions.py & multiGPU-generate_scores_for_unshuffled_insertions.py**: Run the above scripts on multiple GPUs.
- **generate_scores_for_shuffled_insertions.sh & generate_scores_for_unshuffled_insertions.sh**: Automate generating signal strength scores under the slurm system.