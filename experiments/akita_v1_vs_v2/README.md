# AkitaV1 vs. AkitaV2

This directory contains data, matrices, plots, and notebooks related to the AkitaV1 and AkitaV2 prediction analysis.

## Directories:

### 1. data/
This directory contains the following files:
- **v1_v2_sequences.tsv**: Contains all the test sequences overlapping AkitaV1’s and AkitaV2’s test genomic windows.
- **v1_results.tsv, v2_results.tsv**: Include Spearman correlations and MSE between predicted maps and targets.

### 2. matrices/
Contains saved predicted maps as .npz matrices for two genomic windows, for plotting purposes.

### 3. plots/
Includes saved plots in PDF format.

## Notebooks:

- **test_v1_v2_overlap.ipynb**: Generates a TSV table with test genomic windows overlapping both AkitaV1’s and AkitaV2’s test sets.
- **v1.ipynb and v2.ipynb**: Run AkitaV1 and AkitaV2 on test genomic windows.
- **v2_maps.ipynb**: Generates maps for two genomic windows using AkitaV2 for illustration purposes.
- **comparative_analysis.ipynb**: Compares the quality of AkitaV1’s and AkitaV2’s predictions.
- **log-obs-exp_for_targets.ipynb**: Generates log(obs/exp) maps for targets.
- **plotting_maps.ipynb**: Plots maps for two genomic windows to visually compare predictions of AkitaV1 and AkitaV2 with targets.