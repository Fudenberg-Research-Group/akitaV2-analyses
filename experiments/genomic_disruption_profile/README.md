# Genomic Disruption - Profile

This directory contains sliding genomic disruption experiments.

## Directories:

### 1. analysis/
Contains notebooks and plots related to cross-species disruption score analysis.
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing disruption scores, their correlation between human and mouse.

## Files:

- **genomic_disruption_sliding_profile.py**: Script that performs disruption of permutation given a TSV with mouse or human genomic windows, returns disruption scores (SCD) in h5 format.
- **multiGPU_genomic_sliding_profile.py**: Runs the script above on multiple GPUs.
- **genomic_disruption_sliding_profile.sh**: Automates generating scores under the slurm system.