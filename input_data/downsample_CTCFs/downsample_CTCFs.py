#!/usr/bin/env python

import pandas as pd
import numpy as np

import sys
sys.path.insert(0, "/home1/smaruj/akitaX1-analyses/utils/")
from utils import (read_and_average_virtual_exp, read_and_average_genomic_exp)

# stat metric downsampling is based on
stat_of_analysis_interest = "SCD"

# reading data, averaging over targets and backgrounds
insertion_df = read_and_average_virtual_exp("/project/fudenber_735/akitaX1_analyses_data/virtual_insertion_singletons", stat_to_average=stat_of_analysis_interest)

# reading data, averaging over targets
disruption_df = read_and_average_genomic_exp("/project/fudenber_735/akitaX1_analyses_data/genomic_disruption/disruption_by_permutation", stat_to_average=stat_of_analysis_interest)

# ranaming SCD columns
disruption_df = disruption_df.rename(columns={"SCD": "disruption_SCD"})
insertion_df = insertion_df.rename(columns={"SCD": "insertion_SCD"})

df_collected = pd.concat(
    [
        disruption_df["chrom"],
        disruption_df["end"],
        disruption_df["start"],
        disruption_df["strand"],
        disruption_df["disruption_SCD"],
        insertion_df["insertion_SCD"]
    ],
    axis=1,
)

# collecting top1250 sites
top1250 = df_collected.sort_values(by="disruption_SCD", ascending=False)[:1250]
sample250 = df_collected.sort_values(by="disruption_SCD", ascending=False)[1250:].sample(n=250)

subsampled = pd.concat([top1250, sample250])
subsampled = subsampled.reset_index()
subsampled = subsampled.rename(columns={"index": "initial_df_index"})

#saving tsv
subsampled.to_csv("./output/CTCFs_jaspar_filtered_mm10_subsampled.tsv", sep="\t") 
