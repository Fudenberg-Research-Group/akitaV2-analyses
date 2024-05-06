# Virtual Insertion - Flank Length

This directory contains experiments related to virtual insertion analysis with variable flank length.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **analysis_single_insertion.ipynb**: Notebook analyzing insertion of single CTCF sites with variable flank length.
- **single_grid_plotting.ipynb**: Notebook plotting a grid of maps with changing flank length for a randomly chosen CTCF from the class with the highest insertion scores (random_10sites_top_class.tsv).
- **analysis_double_insertion.ipynb**: Notebook analyzing insertion of pairs of CTCF sites with variable flank length.
- **double_grid_plotting.ipynb**: Notebook plotting a grid of maps with changing flank length for CTCF sites with the highest insertion scores (top10_insSCD.tsv).

### 2. input data/
Contains 10 TSV files (one for each background) with a set of mouse CTCF sites.

## Files:

- **generate_insertion_df.py**: Script generating TSV with CTCF sites with insertion parameters.
- **virtual_insertion_with_flanks.py**: Script performing insertion of CTCF sites as specified in the provided TSV table, returns insertion scores in h5 file format.
- **multiGPU_virtual_insertion_with_flanks.py**: Runs the script above on multiple GPUs.
- **virtual_insertion_with_flanks.sh**: Automates generating insertion scores under the Slurm system.