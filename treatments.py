# extract treatment regimens of patients

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
clinical_data = pd.read_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/clinical.tsv', sep='\t')

# keep only case id and treatments.regimen_or_line_of_therapy
treatment_data = clinical_data[['cases.submitter_id', 'treatments.regimen_or_line_of_therapy']]

treatment_data = treatment_data.rename(columns={'cases.submitter_id': 'sample_id', 'treatments.regimen_or_line_of_therapy': 'treatment_regimen'})

# get unique treatment regimens and frequency of each regimen
treatment_counts = treatment_data['treatment_regimen'].value_counts()

treatments = treatment_data['treatment_regimen'].unique()
print(treatment_counts)

def parse_treatment(row):
    txt = row.lower()
    
    # Risk
    if any(k in txt for k in ["low risk", "lr", "sr-low"]):
        risk = "1"
    elif any(k in txt for k in ["average", "sr-avg"]):
        risk = "2"
    elif any(k in txt for k in ["high", "hr", "vhr"]):
        risk = "3"
    else:
        risk = "0"
    
    # Intensity
    if any(k in txt for k in ["intensified", "augmented", "high arm"]):
        intensity = "2"
    else:
        intensity = "1"
    
    # MTX
    if "capizzi" in txt:
        mtx = "1"
    elif "high dose" in txt or "hd-mtx" in txt:
        mtx = "2"
    else:
        mtx = "0"
    
    # Steroid
    if "prednisone" in txt:
        steroid = "1"
    elif "dexamethasone" in txt:
        steroid = "2"
    else:
        steroid = "0"
    
    # Down syndrome
    ds = "1" if "down syndrome" in txt or "ds" in txt else "0"
    
    return pd.Series([risk, intensity, mtx, steroid, ds])

treatment_data[["risk_group", "intensity", "mtx_type", "steroid_type", "down_syndrome"]] = treatment_data["treatment_regimen"].apply(parse_treatment)
print(treatment_data.head())

treatment_data.to_csv('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/treatment_data.csv', index=False)


# histograms for treatment intensity

myclone = 1 # set to 1 to run with myclone data, 0 for pyclone data
if myclone:
    myclone_results_combined_genes = pd.read_csv('myclone-output_combined/myclone_results_combined_genes.tsv', sep='\t')
    # patient_gene_matrix = pd.read_csv('myclone-output_combined/myclone_patient-gene_matrix.tsv', sep='\t')
    fig_nm = 'myclone_top_genes_cellular_prevalence_histograms_treatment.png'
else:
    myclone_results_combined_genes = pd.read_csv('pyclone_output_combined/pyclone_output_combined_genes.tsv', sep='\t')
    # patient_gene_matrix = pd.read_csv('pyclone_output_combined/pyclone_patient-gene_matrix.tsv', sep='\t')
    fig_nm = 'pyclone_top_genes_cellular_prevalence_histograms_intensity.png'

# add treatment to myclone_results_combined_genes
myclone_results_combined_genes = myclone_results_combined_genes.merge(treatment_data, on='sample_id', how='left')

# # # CREATE HISTOGRAM OF CELLULAR FREQUENCY DISTRIBUTION ACROSS TOP 10 MUTATED GENES
# # # i.e. of the patients who have these mutations, 
# # # what is the distribution of cellular prevalence for each gene? 
# # # Are there genes that tend to be more clonal vs subclonal? 
# # # Are there differences in the distribution of cellular prevalence for these top genes between relapse vs non-relapse patients?

# top_genes = ['KRAS'
#     ,'NRAS'
#     ,'CREBBP'
#     ,'PAX5'
#     ,'ETV6'
#     ,'IKZF1'
#     ,'TP53'
#     ,'PTPN11'
#     ,'FLT3'
#     ,'JAK2'
#     ,'NSD2'
#     ,'SETD2'
#     ,'KMT2D'
# ]
# top_genes_df = myclone_results_combined_genes[myclone_results_combined_genes["Hugo_Symbol"].isin(top_genes)]
# # remove clusters with cellular prevalence of <0 (i.e. not actually detected in sample)
# top_genes_df = top_genes_df[top_genes_df["cellular_prevalence"] >= 0]
# # print(top_genes_df[['sample_id', 'treatment_regimen', 'intensity']].head())

# # separate low and high intensity treatment groups
# top_genes_df_low_intensity = top_genes_df[top_genes_df['intensity'] == 1]
# top_genes_df_high_intensity = top_genes_df[top_genes_df['intensity'] == 2]
# top_genes_df_relapse = top_genes_df[top_genes_df['relapse_status']=='Yes']
# top_genes_df_non_relapse = top_genes_df[top_genes_df['relapse_status']=='No']
# # 2 by 5 figure with histograms of cellular frequency distribution for each gene
# fig, axes = plt.subplots(1, 13, figsize=(20, 13))
# for i, gene in enumerate(top_genes):
#     ax = axes[i]
#     sns.histplot(
#         data=top_genes_df_high_intensity[top_genes_df_high_intensity["Hugo_Symbol"] == gene],
#         x="cellular_prevalence",
#         color="red",
#         label="High Intensity",
#         ax=ax,
#         kde=True,
#         stat="density" # use density to normalize for different sample sizes (2x more non-relapse patients than relapse)
#     )
#     sns.histplot(
#         data=top_genes_df_low_intensity[top_genes_df_low_intensity["Hugo_Symbol"] == gene],
#         x="cellular_prevalence",
#         color="blue",
#         label="Low Intensity",
#         ax=ax,
#         kde=True,
#         stat="density"
#     )
#     ax.set_title(gene)
#     ax.legend()
# plt.tight_layout()
# plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/'+fig_nm)
