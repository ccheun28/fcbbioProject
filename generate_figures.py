import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

output_file_path = "/Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv"
#myclone_output_file_path = "/Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv"
#pyclone_output_file_path = "/Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_output_combined/pyclone_output_combined.tsv"
df = pd.read_csv(output_file_path, sep="\t")

clusters_per_sample = (
    df.groupby("sample_id")["cluster_id"]
    .nunique()
    .reset_index(name="num_clusters")
)

#figure: number of clusters per sample

plt.violinplot(
    clusters_per_sample["num_clusters"],
)
plt.ylabel("Number of Clusters Identified per Sample")
plt.xticks([1], ["Samples"])
plt.title("Number of Clusters Identified Per Sample with MyClone")
plt.savefig("figures/myclone_num_clusters_violin_plot.png")
#plt.show()

#figure: number of clusters per sample, relapse vs non relapse
meta = df[["sample_id", "relapse_status"]].drop_duplicates()

plot_df = clusters_per_sample.merge(meta, on="sample_id")
sns.violinplot(
    data=plot_df,
    x="relapse_status",
    y="num_clusters"
)

plt.title("Number of Clusters Identified per Sample with MyClone by Relapse Status")
plt.xlabel("Relapse Status")
plt.ylabel("Number of Clusters")
plt.savefig("figures/myclone_num_clusters_relapse_violin_plot.png")
#plt.show()

#figure: number of clusters per sample vs MRD
mrd = df[["sample_id", "MRD_EOI_Pct"]].drop_duplicates()
plot_df = clusters_per_sample.merge(mrd, on="sample_id")

plt.figure(figsize=(6,5))

plt.scatter(
    plot_df["num_clusters"],
    plot_df["MRD_EOI_Pct"]
)

plt.xlabel("Number of Clusters per Sample")
plt.ylabel("MRD (%)")
plt.title("MRD vs Clonal Complexity (MyClone)")
plt.savefig("figures/myclone_mrd_vs_clusters_scatter.png")
#plt.show()

#figure: Shannon index per sample
shannon_index = df[["sample_id", "shannon_index"]].drop_duplicates()

plot_df = clusters_per_sample.merge(shannon_index, on="sample_id")

sns.violinplot(
    data=plot_df,
    y="shannon_index"
)

plt.title("Clonal Diversity per Sample with MyClone")
plt.xlabel("Samples")
plt.ylabel("Shannon Index")
plt.savefig("figures/myclone_shannon_violin.png")
#plt.show()

#figure: Shannon index per sample relapse vs non-relapse
shannon_df = df[["sample_id", "shannon_index"]].drop_duplicates()
meta = df[["sample_id", "relapse_status"]].drop_duplicates()

plot_df = shannon_df.merge(meta, on="sample_id")

plt.figure(figsize=(6,5))

sns.violinplot(
    data=plot_df,
    x="relapse_status",
    y="shannon_index"
)

plt.title("Clonal Diversity (Shannon Index) by Relapse Status (MyClone)")
plt.xlabel("Relapse Status")
plt.ylabel("Shannon Index")

plt.savefig("figures/myclone_shannon_by_relapse.png")
#plt.show()

#figure: shannon index per sample vs MRD
mrd = df[["sample_id", "MRD_EOI_Pct"]].drop_duplicates()
plot_df = shannon_df.merge(mrd, on="sample_id")

plt.figure(figsize=(6,5))

plt.scatter(
    plot_df["shannon_index"],
    plot_df["MRD_EOI_Pct"]
)

plt.xlabel("Shannon Index")
plt.ylabel("MRD (%)")
plt.title("MRD vs Clonal Diversity (MyClone)")
plt.savefig("figures/myclone_mrd_vs_shannon_scatter.png")
#plt.show()

#figure: Subclonal vs clonal by CCF
df["clone_type"] = df["cellular_prevalence"].apply(
    lambda x: "clonal" if x >= 0.75 else "subclonal"
)

counts = (
    df.groupby(["sample_id", "clone_type"])
    .size()
    .reset_index(name="n_mutations")
)
counts = counts.merge(meta, on="sample_id")

grouped = (
    counts.groupby(["relapse_status", "clone_type"])["n_mutations"]
    .sum()
    .reset_index()
)

pivot = grouped.pivot(
    index="relapse_status",
    columns="clone_type",
    values="n_mutations"
).fillna(0)

# normalize to 100%
pivot = pivot.div(pivot.sum(axis=1), axis=0) * 100

pivot.plot(
    kind="bar",
    stacked=True,
    figsize=(6,5),
    color=["steelblue", "orange"]
)

plt.ylabel("Percentage of Mutations")
plt.title("Clonal vs Subclonal Mutations by Relapse Status (MyClone)")
plt.legend(title="Clone Type")
plt.xticks(rotation=0)

plt.savefig("figures/myclone_clonal_subclonal_relapse_stacked_bar.png")
plt.show()