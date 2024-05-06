# Dot vs. Boundary Experiment

This directory contains the dot vs. boundary experiment.

## Directories:

### 1. analysis/
Contains notebooks and plots related to analysis:
- **plots/**: Directory with saved plots.
- **dots_boundaries_overlap/**:
    - **checking_dots_boundaries_overlap.ipynb**: Notebook comparing the number of CTCF sites overlapping boundaries and dot anchors and both.
    - **dots_boundaries_overlap_HiC.ipynb**: Notebook visualizing boundaries and dot anchors on HiC maps.
- **boundary_windows_analysis.ipynb**: Notebook analyzing predicted dot- and boundary-formation capabilities for boundary-overlapping CTCF sites.
- **dot_anchors_analysis.ipynb**: Notebook analyzing predicted dot- and boundary-formation capabilities for dot anchors-overlapping CTCF sites.
- **picked_for_plotting.tsv**: Table with CTCF sites spread in the range of predicted dot- and boundary-formation capabilities.
- **maps_grid_plotting.ipynb**: Notebook plotting dot and boundary maps for a set of chosen CTCF sites (picked_for_plotting.tsv).

### 2. input data/
- **boundary_windows/**: TSV tables with CTCF sites overlapping boundaries with parameters for boundary and dot scenarios.
- **dot_windows/**: TSV tables with CTCF sites overlapping dot anchors with parameters for boundary and dot scenarios.

## Files:

- **generate_tsv_dot_boundary_scenario.py**: Script generating TSV with CTCF sites to be inserted in dot and boundary scenarios.
- **virtual_symmetric_experiment_dot_vs_boundaries.py**: Script performing insertion of CTCF sites in boundary or dot scenarios, returns boundary or dot strength scores in h5 format.
- **multiGPU_dot_vs_boundaries.py**: Runs the script above on multiple GPUs.
- **generate_boundary_experiment.sh & generate_dot_experiment.sh**: Automating generating boundary and dot scores under the slurm system.