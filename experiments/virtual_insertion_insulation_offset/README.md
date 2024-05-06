# Virtual Insertion - Insulation Offset

This directory contains experiments related to virtual insertion analysis with insulation offset.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **analysis_orientation.ipynb**: Notebook analyzing insulation offset vs orientation asymmetry of CTCF cluster.
- **maps_plotting.ipynb**: Notebook plotting an averaged map for a set of strong CTCF sites (from top10_insSCD.tsv) together with insulation offset profiles.

### 2. input data/
Contains a TSV file with a set of mouse CTCF sites inserted in 6 copies in different orientations.

## Files:

- **generate_insulation_offset_orientation.py**: Script generating TSV with CTCF sites and insertion parameters.
- **virtual_insulation_offset.py**: Script performing insertion of CTCF sites as specified in the provided TSV table, returns insertion scores in h5 file format.
- **multiGPU_virtual_insulation_offset.py**: Runs the script above on multiple GPUs.
- **virtual_insulation_offset.sh**: Automates generating insertion scores under the Slurm system.