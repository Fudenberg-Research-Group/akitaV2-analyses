
## Vary the Number Virtual Insertion Experiment

1. Generate tsv table
   `python generate_number_insertion_df.py --backgrounds-indices "0,1,2,3,4" --output-filename ./input_data/CTCFs_jaspar_filtered_mm10_number_R20_bg04.tsv`

    `python generate_number_insertion_df.py --backgrounds-indices "0" --max-number-sites 30 --output-filename ./input_data/CTCFs_jaspar_filtered_mm10_number_R30_bg0.tsv`

   `python generate_number_insertion_df.py --backgrounds-indices "2" --max-number-sites 40 --output-filename ./input_data/CTCFs_jaspar_filtered_mm10_number_R40_bg2.tsv`

2. Collecting files
    `python collect_jobs_and_clean.py /scratch2/smaruj/vary_number/right40_bg0_m0 -d /home1/smaruj/akitaX1-analyses/experiments/virtual_insertion_number/input_data/CTCFs_jaspar_filtered_mm10_number_R40_bg0.tsv -v -l`

3. Generating predictions
   Change parameters and files in virtual_number_insertion.sh and sbatch it to run the experiement split into jobs.
   