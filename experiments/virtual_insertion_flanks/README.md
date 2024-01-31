
## Single Experiment

1. generating df
   `python generate_insertion_df.py --input-tsv-file /home1/smaruj/akitaX1-analyses/input_data/downsample_CTCFs/output/CTCFs_jaspar_filtered_mm10_subsampled.tsv --backgrounds-indices "0" --output-filename /home1/smaruj/akitaX1-analyses/experiments/virtual_insertion_flanks/input_data/CTCFs_jaspar_filtered_mm10_single_flanks.tsv`

2. prediction -> run virtual_insertion_with_flanks.sh using paramters of choice

## Double Experiment

1. generating df
   `python generate_insertion_df.py --input-tsv-file /home1/smaruj/akitaX1-analyses/input_data/downsample_CTCFs/output/CTCFs_jaspar_filtered_mm10_subsampled.tsv --backgrounds-indices "0" --orientation-string "><" --output-filename /home1/smaruj/akitaX1-analyses/experiments/virtual_insertion_flanks/input_data/CTCFs_jaspar_filtered_mm10_double_convergent_flanks.tsv`

   `python generate_insertion_df.py --input-tsv-file /home1/smaruj/akitaX1-analyses/input_data/downsample_CTCFs/output/CTCFs_jaspar_filtered_mm10_subsampled.tsv --backgrounds-indices "1" --orientation-string "<<" --output-filename /home1/smaruj/akitaX1-analyses/experiments/virtual_insertion_flanks/input_data/CTCFs_jaspar_filtered_mm10_double_left_flanks_bg1.tsv`

2. prediction -> run virtual_insertion_with_flanks.sh using paramters of choice
