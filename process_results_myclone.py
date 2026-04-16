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

### ## ADD GENE NAMES TO MYCLONE RESULTS

input_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv'
maf_input_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/data/raw'
output_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined_genes.tsv'

# map mutation ID to gene name
def map_mutation_to_gene(maf_df):
    maf_df["mutation_id"] = (
        maf_df["Chromosome"].astype(str)
        + ":"
        + maf_df["Start_Position"].astype(str)
        + ":"
        + maf_df["Reference_Allele"]
        + ">"
        + maf_df["Tumor_Seq_Allele2"]
    )
    mutation_to_gene = maf_df[["mutation_id", "Hugo_Symbol"]]
    return mutation_to_gene

# get unique list of patients in myclone_results_combined
myclone_results_combined_df = pd.read_csv(input_path, sep='\t')
myclone_patients = list(set(myclone_results_combined_df["sample_id"]))

# reformat sample sheet
ss = pd.read_csv('gdc_sample_sheet.2026-04-01.tsv', sep='\t') #CHANGED
ss = ss[['File ID', 'File Name', 'Data Type', 'Case ID']].copy()
ss['Case ID'] = ss['Case ID'].str[:13] # only keep first 13 characters to match sample_id format in myclone results
ss = ss.rename(columns={'Case ID': 'sample_id'}) # rename to match sample_id format in myclone results
ss = ss[ss["sample_id"].isin(myclone_patients)] # only keep samples in myclone results
ss_maf = ss.loc[ss['Data Type'] == 'Masked Somatic Mutation'].copy() # only keep MAF files

gene_dfs = [] # list of dataframes to combine later
for case in myclone_patients:
    info = ss_maf.loc[ss_maf['sample_id'] == case].copy()
    maf_path = os.path.join(maf_input_path, info['File ID'].iloc[0], info['File Name'].iloc[0])
    try:
        maf_df = pd.read_csv(maf_path, sep='\t', compression='gzip', comment='#')
    except Exception as e:
        print(f"Error reading MAF file {maf_path}: {e}, skipping sample {case}")
        continue
    
    mutation_to_gene = map_mutation_to_gene(maf_df)

    case_df = myclone_results_combined_df.loc[myclone_results_combined_df['sample_id'] == case].copy()

    case_df = case_df.merge(mutation_to_gene, left_on='mutation_id', right_on='mutation_id', how='left')
    gene_dfs.append(case_df)

# Combine all gene dataframes
myclone_results_combined_genes = pd.concat(gene_dfs, ignore_index=True)

print(myclone_results_combined_genes.shape)
print(myclone_results_combined_genes.head())

myclone_results_combined_genes.to_csv(output_path, sep='\t', index=False)