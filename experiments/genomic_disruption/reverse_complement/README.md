# Disruption by Permutation - Reverse Complement

This directory contains the disruption by permutation experiment with reverse complement predictions.

## Directories:

### 1. analysis/
Contains notebooks and plots related to disruption score analysis.
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing predicted disruption scores for default-strand and reverse complement.

## Files:

- **genomic_disruption_reverse_complement.py**: Script that performs disruption of permutation using the reverse complement DNA strand, given a TSV with mouse CTCF sites, returns disruption scores (SCD) in h5 format.
- **multiGPU_genomic_disruption_reverse_complement.py**: Runs the script above on multiple GPUs.
- **genomic_disruption_reverse_complement.sh**: Automates generating scores under the slurm system.