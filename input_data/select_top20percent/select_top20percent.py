"""
This script processes a TSV file containing strong CTCF sites, filters the top 20% of sites based on their 
insertion SCD scores, and saves the filtered data to a new TSV file.

Usage:
    This script does not take command-line arguments and should be run as is.

The script performs the following steps:
1. Loads a TSV file containing strong CTCF sites data from the specified path.
2. Checks for and removes an unnamed column (typically an index column) if present.
3. Calculates the number of rows corresponding to the top 20% of the dataset.
4. Sorts the data by 'insertion_SCD' in descending order and selects the top 20% of rows.
5. Resets the index of the filtered dataframe to maintain a continuous index.
6. Saves the resulting top 20% of the data to a new TSV file at the specified output path.

Functions and methods used:
    - pd.read_csv: Reads the TSV file into a pandas DataFrame.
    - drop: Removes the unnamed column if present.
    - sort_values: Sorts the dataframe by 'insertion_SCD'.
    - reset_index: Resets the index of the dataframe.
    - to_csv: Saves the filtered dataframe to a TSV file.

Paths:
    - Input path: /home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv
    - Output path: ./output/CTCFs_jaspar_filtered_mm10_top20percent.tsv
"""

import pandas as pd

path = "/home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv"
strong_sites = pd.read_csv(path, sep="\t")

if "Unnamed: 0" in strong_sites.columns:
    strong_sites = strong_sites.drop(columns=["Unnamed: 0"])

len_20percent_data = int(len(strong_sites) * 0.2)

top20percent = strong_sites.sort_values(by="insertion_SCD", ascending=False)[
    :len_20percent_data
]
top20percent = top20percent.reset_index(drop=True)
top20percent.to_csv(
    "./output/CTCFs_jaspar_filtered_mm10_top20percent.tsv",
    sep="\t",
    header=True,
    index=False,
)
