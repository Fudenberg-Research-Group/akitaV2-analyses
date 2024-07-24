"""
This script processes genomic data to analyze the impact of genomic disruptions and virtual insertions on 
SCD (statistic of analysis interest). It reads and averages data over specified targets and backgrounds, 
renames columns for clarity, and collects and subsamples the top significant sites based on disruption scores. 
The final output is a TSV file with the filtered and subsampled data.

Usage:
    This script does not take command-line arguments and should be run as is.

The script performs the following steps:
1. Defines the statistical metric of interest (SCD) for downsampling.
2. Reads and averages virtual insertion data from the specified directory.
3. Reads and averages genomic disruption data from the specified directory.
4. Renames columns in the dataframes for clarity:
   - 'SCD' in the disruption dataframe is renamed to 'disruption_SCD'.
   - 'SCD' in the insertion dataframe is renamed to 'insertion_SCD'.
5. Combines relevant columns from both dataframes into a single dataframe:
   - 'chrom', 'end', 'start', 'strand' from the insertion dataframe.
   - 'disruption_SCD' from the disruption dataframe.
   - 'insertion_SCD' from the insertion dataframe.
6. Sorts the combined dataframe by 'disruption_SCD' in descending order and collects:
   - The top 1250 sites.
   - A random sample of 250 sites from the remaining data.
7. Concatenates the top 1250 sites and the 250 random samples into a single dataframe.
8. Resets the index of the concatenated dataframe and renames the new index column to 'initial_df_index'.
9. Saves the final dataframe to a TSV file at the specified output path.

Functions used:
    - read_and_average_virtual_exp: Reads and averages virtual insertion data.
    - read_and_average_genomic_exp: Reads and averages genomic disruption data.
"""

import pandas as pd
from helper import read_and_average_virtual_exp, read_and_average_genomic_exp

# stat metric downsampling is based on
stat_of_analysis_interest = "SCD"

# reading data, averaging over targets and backgrounds
insertion_df = read_and_average_virtual_exp(
    "/project/fudenber_735/akitaX1_analyses_data/virtual_insertion_singletons",
    stat_to_average=stat_of_analysis_interest,
)

# reading data, averaging over targets
disruption_df = read_and_average_genomic_exp(
    "/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation",
    stat_to_average=stat_of_analysis_interest,
)

# ranaming SCD columns
disruption_df = disruption_df.rename(columns={"SCD": "disruption_SCD"})
insertion_df = insertion_df.rename(columns={"SCD": "insertion_SCD"})

df_collected = pd.concat(
    [
        insertion_df["chrom"],
        insertion_df["end"],
        insertion_df["start"],
        insertion_df["strand"],
        disruption_df["disruption_SCD"],
        insertion_df["insertion_SCD"],
    ],
    axis=1,
)

# collecting top1250 sites
top1250 = df_collected.sort_values(by="disruption_SCD", ascending=False)[:1250]
sample250 = df_collected.sort_values(by="disruption_SCD", ascending=False)[
    1250:
].sample(n=250)

subsampled = pd.concat([top1250, sample250])
subsampled = subsampled.reset_index()
subsampled = subsampled.rename(columns={"index": "initial_df_index"})

# saving tsv
subsampled.to_csv("./output/CTCFs_jaspar_filtered_mm10_strong.tsv", sep="\t")
