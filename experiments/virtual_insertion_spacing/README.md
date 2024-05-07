# Virtual Insertion - Inter-Motif Spacing

This directory contains experiments related to virtual insertion analysis with inter-motif spacing.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing insertion scores vs. inter-motif distance.
- **spacing_grid_plotting.ipynb**: Notebook plotting a grid of maps (for orientation and inter-motif distance) for 10 randomly chosen sites from the top class (random_10sites.tsv).

### 2. input data/
Contains subdirectories:
- **convergent_orientation/**: Contains genomic coordinates of CTCF sites (top 20% of sites having the highest insertion scores in single-site insertion experiment) and with insertion parameters (among others, convergent orientation,inter-motif distance spread logarithmically – the closest sites are to each other, the more often they are tested).
- **divergent_orientation/**: As above, sites in divergent orientation.
- **left_orientation/**: As above, sites in left orientation.
- **right_orientation/**: As above, sites in right orientation.

## Files:

- **generate_spacing_insertion_df.py**: Script generating TSV with CTCF sites and insertion parameters.
- **virtual_insertion_logspacing.py**: Script performing insertion of CTCF sites as specified in the provided TSV table, returns insertion scores in h5 file format.
- **multiGPU_virtual_insertion_logspacing.py**: Runs the script above on multiple GPUs.
- **virtual_insertion_logspacing.sh**: Automates generating insertion scores under the Slurm system.