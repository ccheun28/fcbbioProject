import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import scipy

## CREATE BINARY PATIENT-LEVEL GENE MUTATION MATRIX
# myclone_results_combined_genes = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined_genes.tsv', sep='\t')
# patient_gene_matrix = myclone_results_combined_genes.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
# patient_gene_matrix = (patient_gene_matrix > 0).astype(int) # convert to binary
# print(patient_gene_matrix.shape)
# print(patient_gene_matrix.head())
# patient_gene_matrix.to_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_patient-gene_matrix.tsv', sep='\t', index=True)

## CREATE HEATMAP OF MUTATED GENES ACROSS RELAPSE VS NON-RELAPSE PATIENTS
patient_gene_matrix = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_patient-gene_matrix.tsv', sep='\t', index_col=0)
# relapse_df = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/relapse_dict.csv', sep=',')
# relapse_df = relapse_df.rename(columns={'caseID': 'sample_id'})

# # Merge mutation data
# full_df = relapse_df.merge(
#     patient_gene_matrix,
#     left_on="sample_id",
#     right_index=True,
#     how="left"
# )

# full_df = full_df.sort_values(by='relapse_status')
# # Row colors
# row_colors = full_df["relapse_status"].map({
#     'Yes': "red",
#     'No': "blue"
# })

# only keep relapse/non-relapse patients for heatmap
# full_df = full_df[
#     full_df["relapse_status"] == 'No'
# ]

# create sorted heatmap of all mutated genes across relapse patients
# sort genes by mutation frequency
# full_df = full_df.drop(columns=['relapse_status', 'sample_id'])
# gene_freq = full_df.sum(axis=0).sort_values(ascending=False)
# print(gene_freq.head(20))
# top_n = 75
# top_genes = gene_freq.head(top_n).index
# sorted_matrix = full_df[top_genes]
# print(sorted_matrix.shape)

# g = sns.clustermap(
#     sorted_matrix,
#     row_cluster=False,
#     col_cluster=False,
#     row_colors=row_colors,
#     cmap="Reds",
#     yticklabels=False
# )

# # Force ALL gene labels to show
# g.ax_heatmap.set_xticks(range(len(sorted_matrix.columns)))
# g.ax_heatmap.set_xticklabels(
#     sorted_matrix.columns,
#     rotation=90,
#     fontsize=8
# )
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_patient-gene_matrix_heatmap_all.png')

## CREATE HEATMAP OF MUTATED GENES SORTED BY PATIENT MRD LEVELS
mrd_df = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv', sep='\t')
mrd_df = mrd_df[["sample_id", "MRD_EOI_Pct"]].drop_duplicates(subset="sample_id")

full_df = mrd_df.merge(
    patient_gene_matrix,
    left_on="sample_id",
    right_index=True,
    how="left"
)

full_df = full_df.sort_values(by='MRD_EOI_Pct')

norm = mcolors.Normalize(
    vmin=full_df["MRD_EOI_Pct"].min(),
    vmax=full_df["MRD_EOI_Pct"].max()
)
cmap = plt.cm.coolwarm
row_colors = full_df["MRD_EOI_Pct"].apply(
    lambda x: cmap(norm(x)) if not np.isnan(x) else (1, 1, 1, 1)
)

# sort genes by mutation frequency
full_df = full_df.drop(columns=['MRD_EOI_Pct', 'sample_id'])
gene_freq = full_df.sum(axis=0).sort_values(ascending=False)
print(gene_freq.head(20))
top_n = 75
top_genes = gene_freq.head(top_n).index
sorted_matrix = full_df[top_genes]
print(sorted_matrix.shape)

g = sns.clustermap(
    sorted_matrix,
    row_cluster=False,
    col_cluster=False,
    row_colors=row_colors,
    cmap="Reds",
    yticklabels=False
)

# Force ALL gene labels to show
g.ax_heatmap.set_xticks(range(len(sorted_matrix.columns)))
g.ax_heatmap.set_xticklabels(
    sorted_matrix.columns,
    rotation=90,
    fontsize=8
)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = plt.colorbar(sm, ax=g.ax_heatmap)
cbar.set_label("MRD_EOI_Pct")

plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_patient-gene_matrix_heatmap_mrd.png')

"""
- find top genes --> ccf
- find ccf <0.2 --> subclonal --> genes
- find ccf >0.2 --> clonal --> genes
"""