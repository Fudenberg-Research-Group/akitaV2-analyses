# Flank-Core Compatibility

This directory contains flank-core compatibility analysis.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing core-flank compatibility data.

### 2. input data/
Contains 10 TSV files (one for each background) with a set of 100 strong, 100 medium, and 100 weak CTCF sites (the same set for all backgrounds).

## Files:

- **sites_selection_and_df_generation.py**: Script generating TSV with 100 strong, 100 medium, and 100 weak CTCF sites.
- **virtual_insertion_flank_core_compatibility.py**: Script performing insertion of all combinations of CTCF cores and flanks in the provided TSV table, returns insertion scores of each core-flank pair insertion.
- **multiGPU_virtual_insertion_flank_core_compatibility.py**: Runs the script above on multiple GPUs.
- **virtual_insertion_flank_core_compatibility.sh**: Automates generating insertion scores under the Slurm system.