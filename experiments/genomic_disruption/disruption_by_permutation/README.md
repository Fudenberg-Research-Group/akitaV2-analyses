# Disruption by Permutation

This directory contains the disruption by permutation experiment.

## Directories:

### 1. analysis/
Contains notebooks and plots related to disruption score analysis.
- **plots/**: Directory with saved plots.
- **analysis.ipynb**: Notebook analyzing disruption scores.
- **analysis_disruption_vs_features.ipynb**: Notebook correlating disruption scores with genomic features (CTCF & RAD21 ChIP-seq signal, PWM score, PhyloP score, SMF data).
- **visualize_flanks.ipynb**: Notebook visualizing DNA matrix and creating a logo for sequences with the highest disruption scores.
- **logo_freq_bound.ipynb**: Notebook visualizing DNA matrix and creating a logo for sequences with the highest frequency of the bound state (SMF data).

## Files:

- **genomic_disruption_by_permutation.py**: Script that performs disruption of permutation given a TSV with mouse CTCF sites, returns disruption scores (SCD) in h5 format.
- **multiGPU_genomic_disruption_by_permutation.py**: Runs the script above on multiple GPUs.
- **genomic_disruption_exp.sh**: Automates generating scores under the slurm system.
- **sonmezer_data_genomic_disruption_exp.sh**: Similar to the above, but automates generating scores for the set of CTCF sites tested by Sonmezer_2021.

Sönmezer C, Kleinendorst R, Imanci D, Barzaghi G, Villacorta L, Schübeler D, Benes V, Molina N, Krebs AR. Molecular Co-occupancy Identifies Transcription Factor Binding Cooperativity In Vivo. Mol Cell. 2021 Jan 21;81(2):255-267.e6. doi: 10.1016/j.molcel.2020.11.015. Epub 2020 Dec 7. PMID: 33290745; PMCID: PMC7612519.