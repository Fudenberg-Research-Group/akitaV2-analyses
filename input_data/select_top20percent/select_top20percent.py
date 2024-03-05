
import pandas as pd

path = "/home1/smaruj/akitaX1-analyses/input_data/select_strong_CTCFs/output/CTCFs_jaspar_filtered_mm10_strong.tsv"
strong_sites = pd.read_csv(path, sep="\t")

if "Unnamed: 0" in strong_sites.columns:
    strong_sites = strong_sites.drop(columns=["Unnamed: 0"])

len_20percent_data = int(len(strong_sites) * 0.2)

top20percent = strong_sites.sort_values(by="insertion_SCD", ascending=False)[:len_20percent_data]
top20percent = top20percent.reset_index(drop=True)
top20percent.to_csv("./output/CTCFs_jaspar_filtered_mm10_top20percent.tsv", sep="\t", header=True, index=False)
