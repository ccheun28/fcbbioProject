import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import scipy
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2_contingency



from itertools import combinations
from collections import Counter, defaultdict
import networkx as nx


"""
OVERVIEW OF ALL GENES
"""
## CREATE BINARY PATIENT-LEVEL GENE MUTATION MATRIX
myclone_results_combined_genes = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/pyclone_output_combined/pyclone_output_combined_genes.tsv', sep='\t')
# patient_gene_matrix = myclone_results_combined_genes.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
# patient_gene_matrix = (patient_gene_matrix > 0).astype(int) # convert to binary
# print(patient_gene_matrix.shape)
# print(patient_gene_matrix.head())
# patient_gene_matrix.to_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/pyclone_output_combined/pyclone_patient-gene_matrix.tsv', sep='\t', index=True)

# ## CREATE HEATMAP OF MUTATED GENES ACROSS RELAPSE VS NON-RELAPSE PATIENTS
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
#     full_df["relapse_status"] == 'Yes'
# ]

# create sorted heatmap of all mutated genes across relapse patients
# sort genes by mutation frequency
# full_df = full_df.drop(columns=['relapse_status', 'sample_id'])
# gene_freq = full_df.sum(axis=0).sort_values(ascending=False)
# # print(gene_freq.head(20))
# top_n = 50
# top_genes = gene_freq.head(top_n).index
# print(top_genes)
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
# mrd_df = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/pyclone_output_combined/pyclone_output_combined.tsv', sep='\t')
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

# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_patient-gene_matrix_heatmap_mrd.png')

"""
MORE DETAILED ANALYSIS OF MUTATED GENES IN MYCLONE RESULTS
- find top genes --> ccf
- find ccf <0.2 --> subclonal --> genes
- find ccf >0.2 --> clonal --> genes
- are there genes mutated together/first?
"""

# # CREATE HISTOGRAM OF CELLULAR FREQUENCY DISTRIBUTION ACROSS TOP 10 MUTATED GENES
# # i.e. of the patients who have these mutations, 
# # what is the distribution of cellular prevalence for each gene? 
# # Are there genes that tend to be more clonal vs subclonal? 
# # Are there differences in the distribution of cellular prevalence for these top genes between relapse vs non-relapse patients?

top_genes = ['NRAS'
    ,'PTPN11'
    ,'FLT3'
    ,'JAK2'
    ,'KMT2D'
    ,'NSD2'
]
top_genes_df = myclone_results_combined_genes[myclone_results_combined_genes["Hugo_Symbol"].isin(top_genes)]
# remove clusters with cellular prevalence of <0 (i.e. not actually detected in sample)
top_genes_df = top_genes_df[top_genes_df["cellular_prevalence"] >= 0]

# separate relapse vs non-relapse patients
top_genes_df_relapse = top_genes_df[top_genes_df['relapse_status']=='Yes']
top_genes_df_non_relapse = top_genes_df[top_genes_df['relapse_status']=='No']
# 2 by 5 figure with histograms of cellular frequency distribution for each gene
fig, axes = plt.subplots(1, 6, figsize=(20, 5))
for i, gene in enumerate(top_genes):
    ax = axes[i]
    sns.histplot(
        data=top_genes_df_relapse[top_genes_df_relapse["Hugo_Symbol"] == gene],
        x="cellular_prevalence",
        color="red",
        label="Relapse",
        ax=ax,
        kde=True,
        stat="density" # use density to normalize for different sample sizes (2x more non-relapse patients than relapse)
    )
    sns.histplot(
        data=top_genes_df_non_relapse[top_genes_df_non_relapse["Hugo_Symbol"] == gene],
        x="cellular_prevalence",
        color="blue",
        label="Non-Relapse",
        ax=ax,
        kde=True,
        stat="density"
    )
    ax.set_title(gene)
    ax.legend()
plt.tight_layout()
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_top_genes_cellular_prevalence_histograms_drivers2.png')

"""
Make table of top drivers ccf distribution (clonal, high subclonal, subclonal) in relapse vs non-relapse patients
Proportion of relapse who have KRAS clonal, KRAS subclonal, etc.
"""
top_genes = ['KRAS'
    ,'NRAS'
    ,'CREBBP'
    ,'PAX5'
    ,'ETV6'
    ,'IKZF1'
    ,'TP53'
    ,'PTPN11'
    ,'FLT3'
    ,'JAK2'
    ,'NSD2'
    ,'SETD2'
    ,'KMT2D'
    # ,'CTCF'
    # ,'TRRAP'
    # ,'NCOR2'
    # ,'NT5C2'
    # ,'TBL1XR1'
    # ,'NR3C1'
]

# num_relapse = myclone_results_combined_genes[myclone_results_combined_genes['relapse_status']=='Yes']["sample_id"].nunique()
# num_non_relapse = myclone_results_combined_genes[myclone_results_combined_genes['relapse_status']=='No']["sample_id"].nunique()
# # Categorize each row as clonal (>0.7), high subclonal (0.2-0.7), or subclonal (<0.2)
# myclone_results_combined_genes["clonality_category"] = pd.cut(
#     myclone_results_combined_genes["cellular_prevalence"],
#     bins=[0, 0.2, 0.7, 1],
#     labels=["Subclonal", "High Subclonal", "Clonal"]
# )
# summary_list = []
# for i in top_genes:
#     for category in ["Clonal", "High Subclonal", "Subclonal"]:
#         for relapse_status in ["Yes", "No"]:
#             count = myclone_results_combined_genes[
#                 (myclone_results_combined_genes["Hugo_Symbol"] == i) &
#                 (myclone_results_combined_genes["clonality_category"] == category) &
#                 (myclone_results_combined_genes["relapse_status"] == relapse_status)
#             ]["sample_id"].nunique()
#             total = num_relapse if relapse_status == "Yes" else num_non_relapse
#             proportion = count / total if total > 0 else 0
#             summary_list.append({
#                 "Gene": i,
#                 "Clonality_Category": category,
#                 "Relapse_Status": relapse_status,
#                 "Count": count,
#                 "Proportion": proportion
#             })
# print(summary_list)
# # create grouped bar plot of proportion of patients with each gene mutated in clonal vs subclonal vs high subclonal mutations, separated by relapse vs non-relapse status
# summary_df = pd.DataFrame(summary_list)
# summary_df_pivot = summary_df.pivot_table(
#     index=["Gene", "Relapse_Status"],
#     columns="Clonality_Category",
#     values="Proportion",
#     fill_value=0
# ).reset_index()

# genes = summary_df_pivot["Gene"].unique()
# categories = ["Clonal", "High Subclonal", "Subclonal"]

# x = np.arange(len(genes))  # gene positions
# bar_width = 0.35

# relapse = summary_df_pivot[summary_df_pivot["Relapse_Status"] == "Yes"]
# non_relapse = summary_df_pivot[summary_df_pivot["Relapse_Status"] == "No"]
# relapse = relapse.set_index("Gene").reindex(genes).fillna(0)
# non_relapse = non_relapse.set_index("Gene").reindex(genes).fillna(0)

# plt.figure(figsize=(15, 6))
# fig, ax = plt.subplots(figsize=(15, 4))

# bottom_r = np.zeros(len(genes))
# bottom_nr = np.zeros(len(genes))

# colors = {
#     "Clonal": "#d73027",
#     "High Subclonal": "#fc8d59",
#     "Subclonal": "#91bfdb"
# }

# # relapse bars
# for cat in categories:
#     ax.bar(
#         x - bar_width/2,
#         relapse[cat].values,
#         bar_width,
#         bottom=bottom_r,
#         label=f"{cat} (Relapse)" if x[0] == 0 else "",
#         color=colors[cat]
#     )
#     bottom_r += relapse[cat].values

# # non-relapse bars
# for cat in categories:
#     ax.bar(
#         x + bar_width/2,
#         non_relapse[cat].values,
#         bar_width,
#         bottom=bottom_nr,
#         label=f"{cat} (Non-relapse)" if x[0] == 0 else "",
#         color=colors[cat],
#         alpha=0.6
#     )
#     bottom_nr += non_relapse[cat].values

# ax.set_xticks(x)
# ax.set_xticklabels(genes, rotation=45, ha="right")

# ax.set_ylabel("Proportion of patients")
# ax.set_title("Clonality Distribution by Gene and Relapse Status")

# ax.legend()
# plt.tight_layout()
# # plt.show()
# # plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_top_genes_clonality_distribution.png')

# # Fisher's exact test
# results = []

# for gene in summary_df["Gene"].unique():
    
#     sub = summary_df[summary_df["Gene"] == gene].copy()
    
#     # Collapse into binary clonality
#     sub["Binary_Clonality"] = sub["Clonality_Category"].apply(
#         lambda x: "Clonal" if x == "Clonal" else "Subclonal"
#     )
    
#     # Aggregate counts (IMPORTANT)
#     sub_agg = sub.groupby(
#         ["Relapse_Status", "Binary_Clonality"],
#         as_index=False
#     )["Count"].sum()
    
#     # Build contingency table
#     table = sub_agg.pivot_table(
#         index="Relapse_Status",
#         columns="Binary_Clonality",
#         values="Count",
#         fill_value=0
#     )
    
#     # Ensure consistent structure
#     table = table.reindex(index=["Yes", "No"], fill_value=0)
#     table = table.reindex(columns=["Clonal", "Subclonal"], fill_value=0)

#     # Extract values
#     a = table.loc["Yes", "Clonal"]
#     b = table.loc["Yes", "Subclonal"]
#     c = table.loc["No", "Clonal"]
#     d = table.loc["No", "Subclonal"]
    
#     # Skip if all zeros (no data)
#     if (a + b + c + d) == 0:
#         continue
    
#     # Fisher's exact test
#     odds_ratio, p_value = fisher_exact([[a, b], [c, d]])

#     results.append({
#         "Gene": gene,
#         "Clonal_Relapse": a,
#         "Subclonal_Relapse": b,
#         "Clonal_NonRelapse": c,
#         "Subclonal_NonRelapse": d,
#         "Odds_Ratio": odds_ratio,
#         "p_value": p_value
#     })


# results_df = pd.DataFrame(results)
# results_df["p_adj"] = multipletests(
#     results_df["p_value"],
#     method="fdr_bh"
# )[1]
# # significant = results_df[results_df["p_adj"] < 0.05]
# print(results_df.sort_values("p_value"))



# # FIND TOP GENES IN SUBCLONAL MUTATIONS VS CLONAL MUTATIONS
# # create new patient-gene matrix only for cellular prevalence <0.2 (subclonal)
# subclonal_genes_df = myclone_results_combined_genes[myclone_results_combined_genes["cellular_prevalence"] <= 0.2]
# subclonal_patient_gene_matrix = subclonal_genes_df.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
# subclonal_patient_gene_matrix = (subclonal_patient_gene_matrix > 0).astype(int) # convert to binary
# print(subclonal_patient_gene_matrix.shape)
# print(subclonal_patient_gene_matrix.head())
# subclonal_gene_freq = subclonal_patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
# print(subclonal_gene_freq.head(20))

# # create new patient-gene matrix only for cellular prevalence >0.2 and <0.7 (high subclonal)
# high_subclonal_genes_df = myclone_results_combined_genes[(myclone_results_combined_genes["cellular_prevalence"] > 0.2) & (myclone_results_combined_genes["cellular_prevalence"] <= 0.7)]
# high_subclonal_patient_gene_matrix = high_subclonal_genes_df.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
# high_subclonal_patient_gene_matrix = (high_subclonal_patient_gene_matrix > 0).astype(int)
# print(high_subclonal_patient_gene_matrix.shape)
# print(high_subclonal_patient_gene_matrix.head())
# high_subclonal_gene_freq = high_subclonal_patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
# print(high_subclonal_gene_freq.head(20))

# # create new patient-gene matrix only for cellular prevalence >0.7 (clonal)
# clonal_genes_df = myclone_results_combined_genes[myclone_results_combined_genes["cellular_prevalence"] > 0.7]
# clonal_patient_gene_matrix = clonal_genes_df.groupby(['sample_id', 'Hugo_Symbol']).size().unstack(fill_value=0)
# clonal_patient_gene_matrix = (clonal_patient_gene_matrix > 0).astype(int)
# print(clonal_patient_gene_matrix.shape)
# print(clonal_patient_gene_matrix.head())
# clonal_gene_freq = clonal_patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
# print(clonal_gene_freq.head(20))

# # CREATE HEATMAP OF TOP 20 GENES IN SUBCLONAL MUTATIONS ACROSS RELAPSE VS NON-RELAPSE PATIENTS
# relapse_df_subclonal = relapse_df[relapse_df["sample_id"].isin(subclonal_patient_gene_matrix.index)]
# full_df_subclonal = relapse_df_subclonal.merge(
#     subclonal_patient_gene_matrix,
#     left_on="sample_id",
#     right_index=True,
#     how="left"
# )
# full_df_subclonal = full_df_subclonal.sort_values(by='relapse_status')
# row_colors = full_df_subclonal["relapse_status"].map({
#     'Yes': "red",
#     'No': "blue"
# })
# full_df_subclonal = full_df_subclonal.drop(columns=['relapse_status', 'sample_id'])
# subclonal_gene_freq = full_df_subclonal.sum(axis=0).sort_values(ascending=False)
# top_subclonal_genes = subclonal_gene_freq.head(20).index
# sorted_matrix_subclonal = full_df_subclonal[top_subclonal_genes]
# g_subclonal = sns.clustermap(
#     sorted_matrix_subclonal,
#     row_cluster=False,
#     col_cluster=False,
#     row_colors=row_colors,
#     cmap="Reds",
#     yticklabels=False
# )
# # Force ALL gene labels to show
# g_subclonal.ax_heatmap.set_xticks(range(len(sorted_matrix_subclonal.columns)))
# g_subclonal.ax_heatmap.set_xticklabels(
#     sorted_matrix_subclonal.columns,
#     rotation=90,
#     fontsize=8
# )
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_top_subclonal_genes_heatmap.png')

# # CREATE HEATMAP OF TOP 20 GENES IN HIGH SUBCLONAL MUTATIONS ACROSS RELAPSE VS NON-RELAPSE PATIENTS
# relapse_df_high_subclonal = relapse_df[relapse_df["sample_id"].isin(high_subclonal_patient_gene_matrix.index)]
# full_df_high_subclonal = relapse_df_high_subclonal.merge(
#     high_subclonal_patient_gene_matrix,
#     left_on="sample_id",
#     right_index=True,
#     how="left"
# )
# full_df_high_subclonal = full_df_high_subclonal.sort_values(by='relapse_status')
# row_colors = full_df_high_subclonal["relapse_status"].map({
#     'Yes': "red",
#     'No': "blue"
# })
# full_df_high_subclonal = full_df_high_subclonal.drop(columns=['relapse_status', 'sample_id'])
# high_subclonal_gene_freq = full_df_high_subclonal.sum(axis=0).sort_values(ascending=False)
# top_high_subclonal_genes = high_subclonal_gene_freq.head(20).index
# sorted_matrix_high_subclonal = full_df_high_subclonal[top_high_subclonal_genes]
# g_high_subclonal = sns.clustermap(
#     sorted_matrix_high_subclonal,
#     row_cluster=False,
#     col_cluster=False,
#     row_colors=row_colors,
#     cmap="Reds",
#     yticklabels=False
# )
# # Force ALL gene labels to show
# g_high_subclonal.ax_heatmap.set_xticks(range(len(sorted_matrix_high_subclonal.columns)))
# g_high_subclonal.ax_heatmap.set_xticklabels(
#     sorted_matrix_high_subclonal.columns,
#     rotation=90,
#     fontsize=8
# )
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_top_high_subclonal_genes_heatmap.png')

# # CREATE HEATMAP OF TOP 20 GENES IN CLONAL MUTATIONS ACROSS RELAPSE VS NON-RELAPSE PATIENTS
# relapse_df_clonal = relapse_df[relapse_df["sample_id"].isin(clonal_patient_gene_matrix.index)]
# full_df_clonal = relapse_df_clonal.merge(
#     clonal_patient_gene_matrix,
#     left_on="sample_id",
#     right_index=True,
#     how="left"
# )
# full_df_clonal = full_df_clonal.sort_values(by='relapse_status')
# row_colors = full_df_clonal["relapse_status"].map({
#     'Yes': "red",
#     'No': "blue"
# })
# full_df_clonal = full_df_clonal.drop(columns=['relapse_status', 'sample_id'])
# clonal_gene_freq = full_df_clonal.sum(axis=0).sort_values(ascending=False)
# top_clonal_genes = clonal_gene_freq.head(20).index
# sorted_matrix_clonal = full_df_clonal[top_clonal_genes]
# g_clonal = sns.clustermap(
#     sorted_matrix_clonal,
#     row_cluster=False,
#     col_cluster=False,
#     row_colors=row_colors,
#     cmap="Reds",
#     yticklabels=False
# )
# # Force ALL gene labels to show
# g_clonal.ax_heatmap.set_xticks(range(len(sorted_matrix_clonal.columns)))
# g_clonal.ax_heatmap.set_xticklabels(
#     sorted_matrix_clonal.columns,
#     rotation=90,
#     fontsize=8
# )
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_top_clonal_genes_heatmap.png')

"""
find co-mutated genes in subclonal vs clonal mutations - 20x20 matrices
then co-occurrence networks?
"""

# pair_stats = defaultdict(lambda: {
#     "count": 0,
#     "cellular_prevalence_values": []
# })

# # remove large passenger genes
# remove_genes = [
#     "TTN", "MUC16", "LRP1B", "SYNE1", "OBSCN", "ADGRV1", "ZFHX4", "CSMD3", "DST"
# ]
# remove immunoglobulin genes
# df = myclone_results_combined_genes[~myclone_results_combined_genes["Hugo_Symbol"].str.startswith(("IGH", "IGK", "IGL", "FAT", "DNAH", "MUC"))]
# df = df[~df["Hugo_Symbol"].isin(remove_genes)]

# # keep relapse only
# myclone_results_relapse = df[df['relapse_status']=='No']
# samples = set(myclone_results_relapse["sample_id"])
# for p in samples:
#     df_p = myclone_results_relapse[myclone_results_relapse["sample_id"] == p]
#     for clone_id, group in df_p.groupby("cluster_id"):
#         genes = set(group["Hugo_Symbol"].dropna())
#         clone_ccf = group["cellular_prevalence"].mean()
#         if clone_ccf < 0:
#             continue
#         for g1, g2 in combinations(sorted(genes), 2):
#             pair = (g1, g2)

#             pair_stats[pair]["count"] += 1
#             pair_stats[pair]["cellular_prevalence_values"].append(clone_ccf)

# results = []
# for (g1, g2), stats in pair_stats.items():
#     ccf_array = np.array(stats["cellular_prevalence_values"])

#     results.append({
#         "gene1": g1,
#         "gene2": g2,
#         "count": stats["count"],
#         "mean_ccf": ccf_array.mean(),
#         "median_ccf": np.median(ccf_array),
#         "high_ccf_fraction": np.mean(ccf_array > 0.7),   # clonal fraction
#         "mid_ccf_fraction": np.mean((ccf_array > 0.2) & (ccf_array <= 0.7)), # high subclonal fraction
#         "low_ccf_fraction": np.mean(ccf_array < 0.2)     # subclonal fraction
#     })

# pair_df = pd.DataFrame(results).sort_values("count", ascending=False)
# pair_df = pair_df[ # remove noise by only keeping pairs that co-occur in at least 3 patients and remove pairs of the same gene
#     (pair_df["count"] >= 1) &
#     (pair_df["gene1"] != pair_df["gene2"])
# ]
# print(pair_df.head(20))
# save to csv
# pair_df.to_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/pyclone_output_combined/pyclone_gene_pair_cooccurrence_relapse.tsv', sep='\t', index=False)

# # SCATTERPLOT OF MEAN CCF OF GENE PAIR VS CO-OCCURRENCE COUNT
# plt.figure()
# sns.scatterplot(
#     data=pair_df,
#     x="count",
#     y="mean_ccf"
# )
# plt.xlabel("Co-occurrence frequency")
# plt.ylabel("Mean CCF")
# plt.title("Clonal vs Subclonal Gene Pairs in Relapse Patients")
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_gene_pair_cooccurrence_scatterplot_relapse.png')
# print('scatterplot done')

# # HEATMAP OF TOP 20 GENE PAIRS BY CO-OCCURRENCE COUNT
# top_pairs = pair_df.head(20)
# pivot_count = top_pairs.pivot(index = "gene1", columns = "gene2", values = "count").fillna(0)
# pivot_ccf = top_pairs.pivot(index = "gene1", columns = "gene2", values = "mean_ccf").fillna(0)
# plt.figure(figsize = (15,10))
# sns.heatmap(pivot_count, cmap="Reds")
# sns.heatmap(
#     pivot_count,
#     cmap="Reds",
#     annot=pivot_ccf.round(2),
#     fmt=".2g"
# )
# plt.title("Co-occurrence (color) + CCF (text)")
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_gene_pair_cooccurrence_heatmap_relapse.png')
# print('heatmap done')

# # CREATE CO-OCCURRENCE NETWORK OF TOP GENE PAIRS
# G = nx.Graph()

# # Add edges
# for _, row in pair_df.iterrows():
#     if row["count"] < 3:  # filter noise
#         continue

#     G.add_edge(
#         row["gene1"],
#         row["gene2"],
#         weight=row["count"],
#         ccf=row["mean_ccf"]
#     )

# # Layout
# pos = nx.spring_layout(G, seed=42)

# # Edge properties
# edges = G.edges(data=True)

# weights = [e[2]["weight"] for e in edges]
# ccfs = [e[2]["ccf"] for e in edges]

# # Normalize for plotting
# norm_ccf = [(c - min(ccfs)) / (max(ccfs) - min(ccfs)) for c in ccfs]

# # Draw
# plt.figure(figsize = (15, 15))
# nx.draw_networkx_nodes(G, pos, node_size=200)

# nx.draw_networkx_edges(
#     G, pos,
#     width=[w for w in weights],
#     edge_color=norm_ccf,
#     edge_cmap=plt.cm.coolwarm  # blue → red
# )
# # red edge --> high ccf
# # blue edge --> low ccf
# # thicker --> more frequent cooccurence
# nx.draw_networkx_labels(G, pos)

# plt.title("Gene Co-occurrence Network (CCF encoded)")
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/pyclone_gene_pair_cooccurrence_network_non-relapse_clean.png')