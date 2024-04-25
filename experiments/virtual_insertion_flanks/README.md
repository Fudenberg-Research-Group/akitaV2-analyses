
## Single Experiment

1. generating df
   `python generate_insertion_df.py --input-tsv-file /home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv --backgrounds-indices "0" --output-filename /home1/smaruj/akitaX1-analyses/experiments/virtual_insertion_flanks/input_data/CTCFs_jaspar_filtered_mm10_single_flanks_bg0.tsv`

2. prediction -> run virtual_insertion_with_flanks.sh using paramters of choice

## Double Experiment

1. generating df
   `python generate_insertion_df.py --input-tsv-file /home1/smaruj/akitaX1-analyses/input_data/select_top20percent/output/CTCFs_jaspar_filtered_mm10_top20percent.tsv --orientation-string ">>" --all-permutations --output-filename /home1/smaruj/akitaX1-analyses/experiments/virtual_insertion_flanks/input_data/CTCFs_jaspar_filtered_mm10_top20percent_double.tsv`

2. prediction -> run virtual_insertion_with_flanks.sh using paramters of choice
