"""Perform differential gene expression analysis using pydeseq2."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA
from snakemake.script import snakemake

count_path = Path(snakemake.input["counts"])
samples = pd.read_table(snakemake.input["sample_table"], index_col=0)
conditions = [samples.loc[sample_id, "condition"] for sample_id in samples.index]

if count_path.suffix == ".h5ad":
    ad_counts = ad.read_h5ad(count_path)
    ad_counts.X = ad_counts.X.round().astype(int)
elif count_path.suffix == ".csv":
    df_counts = pd.read_csv(count_path, index_col=0, header=0).sort_index().round().T
    ad_counts = ad.AnnData(X=df_counts)
    ad_counts.obs.index = [index.split("/")[-2] for index in ad_counts.obs.index]
else:
    raise ValueError(f"Unknown file format: {count_path.suffix}")

ad_counts.obs["condition"] = conditions

# plot a PCA to check whether the samples cluster by condition
pca = PCA(n_components=2, svd_solver="full")
pca_counts = pca.fit_transform(ad_counts.X)

if "abundance" in ad_counts.obsm:
    pca_abundances = pca.fit_transform(ad_counts.obsm["abundance"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), dpi=200)
else:
    fig, ax1 = plt.subplots(1, 1, figsize=(5, 5), dpi=200)

sns.scatterplot(x=pca_counts[:, 0], y=pca_counts[:, 1], hue=ad_counts.obs["condition"].values, ax=ax1)
ax1.set_title("PCA of counts")

if "abundance" in ad_counts.obsm:
    sns.scatterplot(x=pca_abundances[:, 0], y=pca_abundances[:, 1], hue=ad_counts.obs["condition"].values, ax=ax2)
    ax2.set_title("PCA of abundances")

fig.tight_layout(pad=5)
fig.savefig(snakemake.output["pca_plot"])

# filter genes with low counts out
if snakemake.params["min_max_counts"] > 0:
    ad_counts = ad_counts[:, ad_counts.X.max(axis=0) >= snakemake.params["min_max_counts"]].copy()

if snakemake.params["min_mean_counts"] > 0:
    ad_counts = ad_counts[:, ad_counts.X.mean(axis=0) >= snakemake.params["min_mean_counts"]].copy()

if snakemake.params["min_sample_counts"]:
    ad_counts = ad_counts[
        :,
        (
            np.sum(ad_counts.X >= snakemake.params["min_sample_counts"]["min_counts"], axis=0)
            >= snakemake.params["min_sample_counts"]["min_samples"]
        ),
    ].copy()

dds = DeseqDataSet(
    adata=ad_counts,
    design_factors="condition",
    refit_cooks=snakemake.params["refit_cooks"],
    inference=DefaultInference(n_cpus=snakemake.threads),
    ref_level=["condition", snakemake.params["untreated_name"]],
    quiet=True,
)

dds.deseq2()

results = DeseqStats(
    dds,
    contrast=(
        "condition",
        snakemake.params["treated_name"],
        snakemake.params["untreated_name"],
    ),
    quiet=True,
)

results.summary()

if snakemake.params["shrink_log2_fold_change"]:
    results.lfc_shrink(coeff=f"condition_{snakemake.params['treated_name']}_vs_{snakemake.params['untreated_name']}")

df_results = results.results_df
df_transcript_to_gene = pd.read_table(
    snakemake.input["transcript_to_gene_map"],
    header=0,
)

# if the transcript to gene map contains a gene name column, replace gene ids with gene names
if "gene_name" in df_transcript_to_gene.columns:
    gene_id_to_gene_name = df_transcript_to_gene.set_index("gene_id")["gene_name"].drop_duplicates().to_dict()
    df_results["gene_name"] = df_results.index.to_series().map(gene_id_to_gene_name).values

# plot a histogram of p-values
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), dpi=200)
sns.histplot(df_results["pvalue"], ax=ax1, bins=20)
ax1.set_title("Histogram of p-values")
sns.histplot(df_results["padj"], ax=ax2, bins=20)
ax2.set_title("Histogram of adjusted p-values")
fig.tight_layout(pad=5)
fig.savefig(snakemake.output["pvalue_histogram"])

df_results.to_csv(
    snakemake.output["results_table"],
    index=True,
)

results.plot_MA(
    log=True,
    save_path=snakemake.output["ma_plot"],
)
