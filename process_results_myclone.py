import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

#for pyclone, can use for myclone with minimal changes

input_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output'
output_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv'

def shannon_index(df):
    cluster_df = df.groupby("cluster_id")["cellular_prevalence"].mean()
    p = cluster_df / cluster_df.sum()
    return -np.sum(p * np.log(p))

pyclone_results_combined = []

for results_folder in os.listdir(input_path):
    path = os.path.join(input_path, results_folder, "myclone_results.tsv")
    df = pd.read_csv(path, sep='\t')
    pyclone_results_combined.append(df)
pyclone_results_combined_df = pd.concat(pyclone_results_combined, ignore_index=True)

relapse_df = pd.read_csv(
    "/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/relapse_dict.csv"
).rename(columns={"caseID": "sample_id"})

mrd_df = (
    pd.read_csv("/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/MP2PRT_ALL_MRD_1510.tsv", sep="\t")
    [["case_submitter_id", "MRD_EOI_Pct"]]
    .rename(columns={"case_submitter_id": "sample_id"})
)

pyclone_results_combined_df["sample_id"] = pyclone_results_combined_df["sample_id"].str.strip()
relapse_df["sample_id"] = relapse_df["sample_id"].str.strip()
mrd_df["sample_id"] = mrd_df["sample_id"].str.strip()

pyclone_results_combined_df = (
    pyclone_results_combined_df
    .merge(mrd_df, on="sample_id", how="left")
    .merge(relapse_df, on="sample_id", how="left")
)

shannon_per_sample = (
    pyclone_results_combined_df
    .groupby("sample_id")
    .apply(shannon_index)
    .reset_index(name="shannon_index")
)

shannon_per_sample["shannon_index"] = shannon_per_sample["shannon_index"].round(4)

pyclone_results_combined_df = pyclone_results_combined_df.merge(
    shannon_per_sample,
    on="sample_id",
    how="left"
)

pyclone_results_combined_df.to_csv(output_path, sep="\t", index=False)

print(pyclone_results_combined_df.shape)
print(pyclone_results_combined_df.head())

