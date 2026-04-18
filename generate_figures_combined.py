import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

myclone_path = "/Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/myclone-output_combined/myclone_results_combined.tsv"
pyclone_path = "/Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_output_combined/pyclone_output_combined.tsv"

os.makedirs("figures", exist_ok=True)

myclone = pd.read_csv(myclone_path, sep="\t")
pyclone = pd.read_csv(pyclone_path, sep="\t")

def get_clusters(df):
    return (
        df.groupby("sample_id")["cluster_id"]
        .nunique()
        .reset_index(name="num_clusters")
    )


my_clusters = get_clusters(myclone)
py_clusters = get_clusters(pyclone)

my_meta = myclone[["sample_id", "relapse_status", "MRD_EOI_Pct", "shannon_index"]].drop_duplicates()
py_meta = pyclone[["sample_id", "relapse_status", "MRD_EOI_Pct", "shannon_index"]].drop_duplicates()


my_clusters = my_clusters.merge(my_meta, on="sample_id")
py_clusters = py_clusters.merge(py_meta, on="sample_id")


#figure: number of clusters per sample
plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
sns.violinplot(y=my_clusters["num_clusters"])
plt.title("MyClone")
plt.ylabel("Num Clusters")

plt.subplot(1,2,2)
sns.violinplot(y=py_clusters["num_clusters"])
plt.title("PyClone")

plt.suptitle("Number of Clusters Identified per Sample")
plt.tight_layout()

plt.savefig("figures/compare_num_clusters_violin.png")
plt.close()


#figure: number of clusters per sample, relapse vs non relapse
plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
sns.violinplot(data=my_clusters, x="relapse_status", y="num_clusters")
plt.title("MyClone")
plt.xlabel("")
plt.ylabel("Num Clusters")

plt.subplot(1,2,2)
sns.violinplot(data=py_clusters, x="relapse_status", y="num_clusters")
plt.title("PyClone")
plt.xlabel("")
plt.ylabel("")

plt.suptitle("Number of Clusters Identified by Relapse Status")
plt.tight_layout()

plt.savefig("figures/compare_clusters_vs_relapse.png")
plt.close()

#figure: number of clusters per sample vs MRD
plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
plt.scatter(my_clusters["num_clusters"], my_clusters["MRD_EOI_Pct"])
plt.title("MyClone")
plt.xlabel("Clusters")
plt.ylabel("MRD")

plt.subplot(1,2,2)
plt.scatter(py_clusters["num_clusters"], py_clusters["MRD_EOI_Pct"])
plt.title("PyClone")
plt.xlabel("Clusters")
plt.ylabel("")

plt.suptitle("MRD vs Number of Clones Identified")
plt.tight_layout()

plt.savefig("figures/compare_mrd_vs_clusters.png")
plt.close()

#figure: Shannon index per sample
plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
sns.violinplot(y=my_clusters["shannon_index"])
plt.title("MyClone")

plt.subplot(1,2,2)
sns.violinplot(y=py_clusters["shannon_index"])
plt.title("PyClone")

plt.suptitle("Clonal Diversity")
plt.tight_layout()

plt.savefig("figures/compare_shannon_violin.png")
plt.close()

#figure: Shannon index per sample relapse vs non-relapse
my_shannon = my_meta[["sample_id", "shannon_index"]]
py_shannon = py_meta[["sample_id", "shannon_index"]]

my_shannon = my_shannon.merge(my_meta[["sample_id", "relapse_status"]], on="sample_id")
py_shannon = py_shannon.merge(py_meta[["sample_id", "relapse_status"]], on="sample_id")

plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
sns.violinplot(data=my_shannon, x="relapse_status", y="shannon_index")
plt.title("MyClone")

plt.subplot(1,2,2)
sns.violinplot(data=py_shannon, x="relapse_status", y="shannon_index")
plt.title("PyClone")

plt.suptitle("Clonal Diversity by Relapse Status")
plt.tight_layout()

plt.savefig("figures/compare_shannon_vs_relapse.png")
plt.close()

#figure: shannon index per sample vs MRD
plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
plt.scatter(my_shannon["shannon_index"], my_clusters["MRD_EOI_Pct"])
plt.title("MyClone")
plt.xlabel("Shannon")
plt.ylabel("MRD")

plt.subplot(1,2,2)
plt.scatter(py_shannon["shannon_index"], py_clusters["MRD_EOI_Pct"])
plt.title("PyClone")
plt.xlabel("Shannon")
plt.ylabel("")

plt.suptitle("MRD vs Shannon Index")
plt.tight_layout()

plt.savefig("figures/compare_mrd_vs_shannon.png")
plt.close()

#figure: Subclonal vs clonal by CCF
def clonal_bar(df, meta, title, ax):
    df = df.copy()
    df["clone_type"] = df["cellular_prevalence"].apply(lambda x: "clonal" if x >= 0.5 else "subclonal")

    counts = (
        df.groupby(["sample_id", "clone_type"])
        .size()
        .reset_index(name="n_mutations")
        .merge(meta, on="sample_id")
    )

    grouped = counts.groupby(["relapse_status", "clone_type"])["n_mutations"].sum().reset_index()

    pivot = grouped.pivot(index="relapse_status", columns="clone_type", values="n_mutations").fillna(0)
    pivot = pivot.div(pivot.sum(axis=1), axis=0) * 100

    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title(title)
    ax.set_ylabel("% Mutations")
    ax.legend().remove()


fig, axes = plt.subplots(1, 2, figsize=(10, 6), constrained_layout=True)

clonal_bar(myclone, my_meta, "MyClone", axes[0])
clonal_bar(pyclone, py_meta, "PyClone", axes[1])

handles, labels = axes[0].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    title="Clone Type",
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    ncol=1
)


fig.suptitle("Percentage of Clonal vs Subclonal Mutations", y=1.05)

plt.savefig(
    "figures/compare_clonal_subclonal.png",
    bbox_inches="tight",
    dpi=300
)

plt.close()