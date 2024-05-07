# Virtual Insertion - Pairwise Mutagenesis

This directory contains experiments related to virtual insertion analysis with pairwise mutagenesis.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing insertion scores in the pairwise mutagenesis experiment.

### 2. input data/
Contains one TSV file per background, with genomic coordinates of CTCF sites (100 CTCFs with the highest insertion score in the single-site insertion experiment) with insertion parameters (among others, pairs of positions to be mutated).

## Files:

- **generate_pairwise_mutagenesis_df.py**: Script generating TSV with CTCF sites and insertion parameters.
- **virtual_insertion_pairwise_mutagenesis.py**: Script performing insertion of CTCF sites as specified in the provided TSV table, returns insertion scores in h5 file format.
- **multiGPU_virtual_insertion_pairwise_mutagenesis.py**: Runs the script above on multiple GPUs.
- **virtual_insertion_pairwise_mutagenesis.sh**: Automates generating insertion scores under the Slurm system.