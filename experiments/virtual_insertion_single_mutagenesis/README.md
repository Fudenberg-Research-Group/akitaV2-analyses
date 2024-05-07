# Virtual Insertion - Single Mutagenesis

This directory contains experiments related to virtual insertion analysis with single mutagenesis.

## Directories:

### 1. analysis/
Contains notebooks and saved data related to analysis:
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing insertion scores in the single mutagenesis experiment.
- **single_mutation_SCD.npz**: Saved insertion scores for original nucleotides at each position, and average insertion scores for substituted nucleotides at each position.

### 2. input data/
Contains a TSV file with genomic coordinates of CTCF sites (100 strong CTCF sites = with the highest insertion scores in single-insertion experiment) with insertion parameters (e.g. positions to be mutated, substitutions).

## Files:

- **generate_single_mutagenesis_df.py**: Script generating TSV with CTCF sites and insertion parameters.
- **virtual_insertion_single_mutagenesis.py**: Script performing insertion of CTCF sites as specified in the provided TSV table, returns insertion scores in h5 file format.
- **multiGPU_virtual_insertion_single_mutagenesis.py**: Runs the script above on multiple GPUs.
- **virtual_insertion_single_mutagenesis.sh**: Automates generating insertion scores under the Slurm system.