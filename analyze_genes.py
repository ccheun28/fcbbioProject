import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import scipy

"""
OVERVIEW OF ALL GENES
"""
## CREATE BINARY PATIENT-LEVEL GENE MUTATION MATRIX
myclone_results_combined_genes = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined_genes.tsv', sep='\t')
# patient_gene_matrix = myclone_results_combined_genes.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
# patient_gene_matrix = (patient_gene_matrix > 0).astype(int) # convert to binary
# print(patient_gene_matrix.shape)
# print(patient_gene_matrix.head())
# patient_gene_matrix.to_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_patient-gene_matrix.tsv', sep='\t', index=True)

## CREATE HEATMAP OF MUTATED GENES ACROSS RELAPSE VS NON-RELAPSE PATIENTS
patient_gene_matrix = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_patient-gene_matrix.tsv', sep='\t', index_col=0)
relapse_df = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/relapse_dict.csv', sep=',')
relapse_df = relapse_df.rename(columns={'caseID': 'sample_id'})

# Merge mutation data
full_df = relapse_df.merge(
    patient_gene_matrix,
    left_on="sample_id",
    right_index=True,
    how="left"
)

full_df = full_df.sort_values(by='relapse_status')
# Row colors
row_colors = full_df["relapse_status"].map({
    'Yes': "red",
    'No': "blue"
})

# only keep relapse/non-relapse patients for heatmap
# full_df = full_df[
#     full_df["relapse_status"] == 'No'
# ]

# create sorted heatmap of all mutated genes across relapse patients
# sort genes by mutation frequency
full_df = full_df.drop(columns=['relapse_status', 'sample_id'])
gene_freq = full_df.sum(axis=0).sort_values(ascending=False)
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
# mrd_df = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv', sep='\t')
# mrd_df = mrd_df[["sample_id", "MRD_EOI_Pct"]].drop_duplicates(subset="sample_id")

# full_df = mrd_df.merge(
#     patient_gene_matrix,
#     left_on="sample_id",
#     right_index=True,
#     how="left"
# )

# full_df = full_df.sort_values(by='MRD_EOI_Pct')

# norm = mcolors.Normalize(
#     vmin=full_df["MRD_EOI_Pct"].min(),
#     vmax=full_df["MRD_EOI_Pct"].max()
# )
# cmap = plt.cm.coolwarm
# row_colors = full_df["MRD_EOI_Pct"].apply(
#     lambda x: cmap(norm(x)) if not np.isnan(x) else (1, 1, 1, 1)
# )

# # sort genes by mutation frequency
# full_df = full_df.drop(columns=['MRD_EOI_Pct', 'sample_id'])
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

# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])

# cbar = plt.colorbar(sm, ax=g.ax_heatmap)
# cbar.set_label("MRD_EOI_Pct")

# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_patient-gene_matrix_heatmap_mrd.png')

"""
MORE DETAILED ANALYSIS OF MUTATED GENES IN MYCLONE RESULTS
- find top genes --> ccf
- find ccf <0.2 --> subclonal --> genes
- find ccf >0.2 --> clonal --> genes
- are there genes mutated together/first?
"""

# CREATE HISTOGRAM OF CELLULAR FREQUENCY DISTRIBUTION ACROSS TOP 10 MUTATED GENES
# i.e. of the patients who have these mutations, 
# what is the distribution of cellular prevalence for each gene? 
# Are there genes that tend to be more clonal vs subclonal? 
# Are there differences in the distribution of cellular prevalence for these top genes between relapse vs non-relapse patients?
# top_genes = gene_freq.head(20).index
# top_genes_df = myclone_results_combined_genes[myclone_results_combined_genes["Hugo_Symbol"].isin(top_genes)]
# # remove clusters with cellular prevalence of <0 (i.e. not actually detected in sample)
# top_genes_df = top_genes_df[top_genes_df["cellular_prevalence"] >= 0]

# # separate relapse vs non-relapse patients
# top_genes_df_relapse = top_genes_df[top_genes_df['relapse_status']=='Yes']
# top_genes_df_non_relapse = top_genes_df[top_genes_df['relapse_status']=='No']
# # 2 by 5 figure with histograms of cellular frequency distribution for each gene
# fig, axes = plt.subplots(4, 5, figsize=(20, 20))
# for i, gene in enumerate(top_genes):
#     ax = axes[i // 5, i % 5]
#     sns.histplot(
#         data=top_genes_df_relapse[top_genes_df_relapse["Hugo_Symbol"] == gene],
#         x="cellular_prevalence",
#         color="red",
#         label="Relapse",
#         ax=ax,
#         kde=True,
#         stat="density" # use density to normalize for different sample sizes (2x more non-relapse patients than relapse)
#     )
#     sns.histplot(
#         data=top_genes_df_non_relapse[top_genes_df_non_relapse["Hugo_Symbol"] == gene],
#         x="cellular_prevalence",
#         color="blue",
#         label="Non-Relapse",
#         ax=ax,
#         kde=True,
#         stat="density"
#     )
#     ax.set_title(gene)
#     ax.legend()
# plt.tight_layout()
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_top_genes_cellular_prevalence_histograms.png')


# FIND TOP GENES IN SUBCLONAL MUTATIONS VS CLONAL MUTATIONS
# create new patient-gene matrix only for cellular prevalence <0.2 (subclonal)
subclonal_genes_df = myclone_results_combined_genes[myclone_results_combined_genes["cellular_prevalence"] <= 0.2]
subclonal_patient_gene_matrix = subclonal_genes_df.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
subclonal_patient_gene_matrix = (subclonal_patient_gene_matrix > 0).astype(int) # convert to binary
print(subclonal_patient_gene_matrix.shape)
print(subclonal_patient_gene_matrix.head())
subclonal_gene_freq = subclonal_patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
print(subclonal_gene_freq.head(20))

# create new patient-gene matrix only for cellular prevalence >0.2 and <0.7 (high subclonal)
high_subclonal_genes_df = myclone_results_combined_genes[(myclone_results_combined_genes["cellular_prevalence"] > 0.2) & (myclone_results_combined_genes["cellular_prevalence"] <= 0.7)]
high_subclonal_patient_gene_matrix = high_subclonal_genes_df.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
high_subclonal_patient_gene_matrix = (high_subclonal_patient_gene_matrix > 0).astype(int)
print(high_subclonal_patient_gene_matrix.shape)
print(high_subclonal_patient_gene_matrix.head())
high_subclonal_gene_freq = high_subclonal_patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
print(high_subclonal_gene_freq.head(20))

# create new patient-gene matrix only for cellular prevalence >0.7 (clonal)
clonal_genes_df = myclone_results_combined_genes[myclone_results_combined_genes["cellular_prevalence"] > 0.7]
clonal_patient_gene_matrix = clonal_genes_df.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
clonal_patient_gene_matrix = (clonal_patient_gene_matrix > 0).astype(int)
print(clonal_patient_gene_matrix.shape)
print(clonal_patient_gene_matrix.head())
clonal_gene_freq = clonal_patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
print(clonal_gene_freq.head(20))

# CREATE HEATMAP OF TOP 20 GENES IN SUBCLONAL MUTATIONS ACROSS RELAPSE VS NON-RELAPSE PATIENTS
relapse_df_subclonal = relapse_df[relapse_df["sample_id"].isin(subclonal_patient_gene_matrix.index)]
full_df_subclonal = relapse_df_subclonal.merge(
    subclonal_patient_gene_matrix,
    left_on="sample_id",
    right_index=True,
    how="left"
)
full_df_subclonal = full_df_subclonal.sort_values(by='relapse_status')
row_colors = full_df_subclonal["relapse_status"].map({
    'Yes': "red",
    'No': "blue"
})
full_df_subclonal = full_df_subclonal.drop(columns=['relapse_status', 'sample_id'])
subclonal_gene_freq = full_df_subclonal.sum(axis=0).sort_values(ascending=False)
top_subclonal_genes = subclonal_gene_freq.head(20).index
sorted_matrix_subclonal = full_df_subclonal[top_subclonal_genes]
g_subclonal = sns.clustermap(
    sorted_matrix_subclonal,
    row_cluster=False,
    col_cluster=False,
    row_colors=row_colors,
    cmap="Reds",
    yticklabels=False
)
# Force ALL gene labels to show
g_subclonal.ax_heatmap.set_xticks(range(len(sorted_matrix_subclonal.columns)))
g_subclonal.ax_heatmap.set_xticklabels(
    sorted_matrix_subclonal.columns,
    rotation=90,
    fontsize=8
)
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_top_subclonal_genes_heatmap.png')

# CREATE HEATMAP OF TOP 20 GENES IN HIGH SUBCLONAL MUTATIONS ACROSS RELAPSE VS NON-RELAPSE PATIENTS
relapse_df_high_subclonal = relapse_df[relapse_df["sample_id"].isin(high_subclonal_patient_gene_matrix.index)]
full_df_high_subclonal = relapse_df_high_subclonal.merge(
    high_subclonal_patient_gene_matrix,
    left_on="sample_id",
    right_index=True,
    how="left"
)
full_df_high_subclonal = full_df_high_subclonal.sort_values(by='relapse_status')
row_colors = full_df_high_subclonal["relapse_status"].map({
    'Yes': "red",
    'No': "blue"
})
full_df_high_subclonal = full_df_high_subclonal.drop(columns=['relapse_status', 'sample_id'])
high_subclonal_gene_freq = full_df_high_subclonal.sum(axis=0).sort_values(ascending=False)
top_high_subclonal_genes = high_subclonal_gene_freq.head(20).index
sorted_matrix_high_subclonal = full_df_high_subclonal[top_high_subclonal_genes]
g_high_subclonal = sns.clustermap(
    sorted_matrix_high_subclonal,
    row_cluster=False,
    col_cluster=False,
    row_colors=row_colors,
    cmap="Reds",
    yticklabels=False
)
# Force ALL gene labels to show
g_high_subclonal.ax_heatmap.set_xticks(range(len(sorted_matrix_high_subclonal.columns)))
g_high_subclonal.ax_heatmap.set_xticklabels(
    sorted_matrix_high_subclonal.columns,
    rotation=90,
    fontsize=8
)
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_top_high_subclonal_genes_heatmap.png')

# CREATE HEATMAP OF TOP 20 GENES IN CLONAL MUTATIONS ACROSS RELAPSE VS NON-RELAPSE PATIENTS
relapse_df_clonal = relapse_df[relapse_df["sample_id"].isin(clonal_patient_gene_matrix.index)]
full_df_clonal = relapse_df_clonal.merge(
    clonal_patient_gene_matrix,
    left_on="sample_id",
    right_index=True,
    how="left"
)
full_df_clonal = full_df_clonal.sort_values(by='relapse_status')
row_colors = full_df_clonal["relapse_status"].map({
    'Yes': "red",
    'No': "blue"
})
full_df_clonal = full_df_clonal.drop(columns=['relapse_status', 'sample_id'])
clonal_gene_freq = full_df_clonal.sum(axis=0).sort_values(ascending=False)
top_clonal_genes = clonal_gene_freq.head(20).index
sorted_matrix_clonal = full_df_clonal[top_clonal_genes]
g_clonal = sns.clustermap(
    sorted_matrix_clonal,
    row_cluster=False,
    col_cluster=False,
    row_colors=row_colors,
    cmap="Reds",
    yticklabels=False
)
# Force ALL gene labels to show
g_clonal.ax_heatmap.set_xticks(range(len(sorted_matrix_clonal.columns)))
g_clonal.ax_heatmap.set_xticklabels(
    sorted_matrix_clonal.columns,
    rotation=90,
    fontsize=8
)
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_top_clonal_genes_heatmap.png')

"""
find co-mutated genes in subclonal vs clonal mutations - 20x20 matrices
then co-occurrence networks?
"""
