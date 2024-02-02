
## Core vs Flank Compatibility Experiment

1. Sample sites and created input df
   Tsv files for all backgrounds have to be generated at once, since the CTCF sites are chosen randomly. 
   `python sites_selection_and_df_generation.py --orientation-string ">" --output-filename CTCFs_jaspar_filtered_mm10_sampled_high_medium_low.tsv`

2. Prediction generation
   After updating parameters, run `virtual_insertion_flank_core_compatibility.sh`
   