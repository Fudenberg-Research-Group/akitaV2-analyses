# Genomic Disruption - Shifted Permutations

This directory contains experiments related to shifted permutations analysis.

## Directories:

### 1. analysis/
Contains notebooks and plots related to disruption score analysis.
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing predicted disruption scores for centered and shifted permutations.

## Files:

- **genomic_shifted_permutation.py**: Script that performs disruption by permutation centering the window on the permuted region, or shifted up to +/-10kb, returns disruption scores (SCD) in h5 format.
- **multiGPU_genomic_shifted_permutation.py**: Runs the script above on multiple GPUs.
- **genomic_shifted_permutation.sh**: Automates generating scores under the slurm system.