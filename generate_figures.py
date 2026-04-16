import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

output_file_path = "/Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_output_combined/pyclone_output_combined.tsv"

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
plt.title("Number of Clusters Identified Per Sample with PyClone")
plt.savefig("figures/pyclone_num_clusters_violin_plot.png")
plt.show()

#figure: number of clusters per sample, relapse vs non relapse
meta = df[["sample_id", "relapse_status"]].drop_duplicates()

plot_df = clusters_per_sample.merge(meta, on="sample_id")
sns.violinplot(
    data=plot_df,
    x="relapse_status",
    y="num_clusters"
)

plt.title("Number of Clusters Identified per Sample with PyClone by Relapse Status")
plt.xlabel("Relapse Status")
plt.ylabel("Number of Clusters")
plt.savefig("figures/pyclone_num_clusters_relapse_violin_plot.png")
plt.show()