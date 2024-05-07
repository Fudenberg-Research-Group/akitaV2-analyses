# Virtual Insertion - Number of Inserts

This directory contains experiments related to virtual insertion analysis with varying numbers of inserts.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing insertion score vs. number of inserts.

### 2. input data/
Contains one TSV file per background, with genomic coordinates of CTCF sites with insertion parameters (up to 40 inserts into background sequence).

## Files:

- **generate_number_insertion_df.py**: Script generating TSV with CTCF sites and insertion parameters.
- **virtual_number_insertion.py**: Script performing insertion of CTCF sites as specified in the provided TSV table, returns insertion scores in h5 file format.
- **multiGPU_virtual_number_insertion.py**: Runs the script above on multiple GPUs.
- **virtual_number_insertion.sh**: Automates generating insertion scores under the Slurm system.   