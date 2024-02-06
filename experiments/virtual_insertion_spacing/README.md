
# Spacing Virtual Insertion Experiment

1. Generate input df with experiment specification
   `python generate_spacing_insertion_df.py --orientation-string ">>" --space-range 1,500000 --num_log-intervals 400 --backgrounds-indices 0 --output-filename ./input_data/CTCFs_jaspar_filtered_mm10_subsampled_logspacing_bg0_R.tsv`

   `python generate_spacing_insertion_df.py --orientation-string "<<" --space-range 1,500000 --num_log-intervals 400 --backgrounds-indices 0 --output-filename ./input_data/CTCFs_jaspar_filtered_mm10_subsampled_logspacing_bg0_L.tsv`

   `python generate_spacing_insertion_df.py --orientation-string "<>" --space-range 1,500000 --num_log-intervals 400 --backgrounds-indices 0 --output-filename ./input_data/CTCFs_jaspar_filtered_mm10_subsampled_logspacing_bg0_D.tsv`

   `python generate_spacing_insertion_df.py --orientation-string "><" --space-range 1,500000 --num_log-intervals 400 --backgrounds-indices 0 --output-filename ./input_data/CTCFs_jaspar_filtered_mm10_subsampled_logspacing_bg0_C.tsv`

2. Generating Predictions with Akita
   Run virtual_insertion_logspacing.sh after changing parameters and input tsv table. 
