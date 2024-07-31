# Input Data

This directory contains input data used in various experiments.

## Directories:

### 1. preprocess_boundary_CTCFs/
Contains a python script and an output TSV with filtered 7,560 mouse (mm10) CTCF sites overlapping TAD boundaries (CTCF sites don’t overlap each other neither repeatable elements).

### 2. preprocess_dot_CTCFs/
Contains a python script and an output TSV with filtered 36,948 mouse (mm10) CTCF sites overlapping dot anchors, but non-overlapping TAD-boundaries (additionally, CTCF sites don’t overlap each other neither repeatable elements).

### 3. select_strong_CTCFs/
Contains a python script and an output TSV with filtered 1,500 mouse (mm10) CTCF sites overlapping TAD boundaries: 1,250 with the highest disruption scores in the disruption by permutation experiment and 250 sites selected randomly from the remaining pool.

### 4. select_top20percent/
Contains a python script and an output TSV with filtered top20% mouse (mm10) CTCF sites overlapping TAD boundaries and yielding the highest insertion scores in the single-site insertion experiment.

### 5. sonmezer2021_SMF_CTCF_binding_data/
Contains R and Jupyter notebooks converting data provided in the Sonmezer2021 paper into a single TSV table.

Sönmezer C, Kleinendorst R, Imanci D, Barzaghi G, Villacorta L, Schübeler D, Benes V, Molina N, Krebs AR. Molecular Co-occupancy Identifies Transcription Factor Binding Cooperativity In Vivo. Mol Cell. 2021 Jan 21;81(2):255-267.e6. doi: 10.1016/j.molcel.2020.11.015. Epub 2020 Dec 7. PMID: 33290745; PMCID: PMC7612519.