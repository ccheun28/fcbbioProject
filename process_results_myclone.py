#import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

#for pyclone, can use for myclone with minimal changes

input_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output'
output_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv'

def shannon_index(df):
    df = df[df["cellular_prevalence"] > 0]

    cluster_df = df.groupby("cluster_id")["cellular_prevalence"].mean()

    if cluster_df.empty or cluster_df.sum() == 0:
        return np.nan

    p = cluster_df / cluster_df.sum()
    return -np.sum(p * np.log(p))

myclone_results_combined = [] # list of dataframes to combine later

for results_folder in os.listdir(input_path):
    path = os.path.join(input_path, results_folder, "myclone_result.tsv")
    try:
        df = pd.read_csv(path, sep='\t')
    except FileNotFoundError:
        #print(f"Warning: {path} not found. Skipping this sample.")
        continue
    except NotADirectoryError:
        #print(f"Error reading {path}: {e}. Skipping this sample.")
        continue
    
    df["sample_id"] = results_folder[-13:] # add sample ID column based on folder name
    df = df.rename(columns={"myclone_id": "cluster_id"}) # rename to match pyclone output format
    df = df.rename(columns={"myclone_ccf": "cellular_prevalence"}) # rename to match pyclone output format
    myclone_results_combined.append(df)
myclone_results_combined_df = pd.concat(myclone_results_combined, ignore_index=True)

## add clinical info
relapse_df = pd.read_csv(
    "/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/relapse_dict.csv"
).rename(columns={"caseID": "sample_id"})

mrd_df = (
    pd.read_csv("/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/MP2PRT_ALL_MRD_1510.tsv", sep="\t")
    [["case_submitter_id", "MRD_EOI_Pct"]]
    .rename(columns={"case_submitter_id": "sample_id"})
)

myclone_results_combined_df["sample_id"] = myclone_results_combined_df["sample_id"].str.strip()
relapse_df["sample_id"] = relapse_df["sample_id"].str.strip()
mrd_df["sample_id"] = mrd_df["sample_id"].str.strip()

myclone_results_combined_df = (
    myclone_results_combined_df
    .merge(mrd_df, on="sample_id", how="left")
    .merge(relapse_df, on="sample_id", how="left")
)

## calculate Shannon index
shannon_per_sample = (
    myclone_results_combined_df
    .groupby("sample_id")
    .apply(shannon_index)
    .reset_index(name="shannon_index")
)

shannon_per_sample["shannon_index"] = shannon_per_sample["shannon_index"].round(4)

myclone_results_combined_df = myclone_results_combined_df.merge(
    shannon_per_sample,
    on="sample_id",
    how="left"
)

myclone_results_combined_df.to_csv(output_path, sep="\t", index=False)

print(myclone_results_combined_df.shape)
print(myclone_results_combined_df.head())

