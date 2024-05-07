# Single-Site Virtual Insertion

This directory contains experiments related to single-site insertion analysis.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing insertion scores in the single-site insertion experiment.
- **insertion_disruption_correlation.ipynb**: Notebook correlating insertion score in this experiment with disruption score from the disruption by permutation experiment.
- **insertion_disruption_correlation-colored_by_number.ipynb**: Copy of the above notebook, with CTCF sites colored by the number of sites in the +-10kb window.
- **visualize_flanks.ipynb**: Notebook analyzing sequence preferences in the flanking regions.
- **downstream_meme.txt & upstream_meme.txt**: Sequence preferences in the downstream and upstream flanks in the MEME format.

### 2. input data/
Contains a TSV file with genomic coordinates of CTCF sites (n=7,560) with insertion parameters (among others, flanking regions length = 30bp).

## Files:

- **generate_single_insertion_df.py**: Script generating TSV with CTCF sites and insertion parameters.
- **virtual_single_insertion.py**: Script performing insertion of CTCF sites as specified in the provided TSV table, returns insertion scores in h5 file format.
- **multiGPU_virtual_single_insertion.py**: Runs the script above on multiple GPUs.
- **virtual_single_insertion.sh**: Automates generating insertion scores under the Slurm system.