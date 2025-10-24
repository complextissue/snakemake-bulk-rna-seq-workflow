"""Perform gene set enrichment analysis base: volcano plot, transcription factor and pathway enrichment analysis."""

from logging import warning
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from snakemake.script import snakemake

df_results = pd.read_csv(snakemake.input["results_table"], index_col=0, header=0)

# replace gene ids with gene names
if "gene_name" in df_results.columns:
    print("Found gene_name column in df_results")
    print(f"Number of unique gene_ids: {df_results.index.nunique()}")
    print(f"Number of unique gene_names: {df_results['gene_name'].nunique()}")
    print(f"Sample gene_ids: {list(df_results.index[:10])}")
    print(f"Sample gene_names: {list(df_results['gene_name'].head(10))}")

    if df_results["gene_name"].nunique() != df_results.index.nunique():
        warning(
            "Missing or duplicated gene names, setting to gene_ids:\n"
            f"{df_results[df_results.duplicated('gene_name')].head(5)}"
        )
        df_results.loc[df_results.duplicated("gene_name"), "gene_name"] = df_results[
            df_results.duplicated("gene_name")
        ].index

    df_results = df_results.set_index("gene_name")
    print(
        f"After setting gene_name as index, sample indices: {list(df_results.index[:10])}"
    )
else:
    print("No gene_name column found, using existing index")
    print(f"Sample gene indices: {list(df_results.index[:10])}")

# save a volcano plot of the differentially expressed genes
fig = dc.pl.volcano(
    data=df_results,
    x="log2FoldChange",
    y="padj",
    thr_stat=snakemake.params["log2_fold_change_threshold"],
    thr_sign=snakemake.params["significance_threshold"],
    top=snakemake.params["top_genes"],
    figsize=(10, 7.5),
    return_fig=True,
)
count_up = np.sum(
    (df_results["log2FoldChange"] > snakemake.params["log2_fold_change_threshold"])
    & (df_results["padj"] < snakemake.params["significance_threshold"])
)
count_down = np.sum(
    (df_results["log2FoldChange"] < -snakemake.params["log2_fold_change_threshold"])
    & (df_results["padj"] < snakemake.params["significance_threshold"])
)
fig.suptitle(
    f"Differentially expressed genes ({count_up} up, {count_down} down)", fontsize=12
)
fig.savefig(snakemake.output["volcano_plot"], dpi=500, bbox_inches="tight")
plt.close(fig)

# perform transcription factor enrichment analysis
collectri = pd.read_csv(snakemake.input["collectri"], index_col=0)

treated_vs_untreated_identifier = f"{snakemake.params['treated_name']}.vs.{snakemake.params['untreated_name']}"

mat = df_results[["stat"]].T.rename(
    index={"stat": treated_vs_untreated_identifier},
)

tf_acts, tf_pvals = dc.mt.ulm(
    data=mat,
    net=collectri,
    verbose=True,
    tmin=10,
)

fig = dc.pl.barplot(
    data=tf_acts,
    name=treated_vs_untreated_identifier,
    top=snakemake.params["top_transcription_factors"],
    vertical=False,
    figsize=(10, 5),
    return_fig=True,
)
fig.tight_layout(pad=1)
transcription_factors_barplot_path = Path(
    snakemake.output["transcription_factors_barplot"]
)
fig.savefig(transcription_factors_barplot_path, dpi=500, bbox_inches="tight")
plt.close(fig)
tf_acts.to_csv(snakemake.output["transcription_factors_table"])

# Plot the log2 fold change of the target genes of the top transcription factors
tf_acts_ranked = tf_acts.T.copy()
tf_acts_ranked["activity_absolute"] = np.abs(
    tf_acts_ranked[treated_vs_untreated_identifier]
)
top_tfs = (
    tf_acts_ranked.sort_values("activity_absolute", ascending=False)
    .head(snakemake.params["top_transcription_factors"])
    .index
)

# Get the output directory for TF-specific plots
tf_plots_dir = transcription_factors_barplot_path.parent

for top_tf in top_tfs:
    # Volcano plot with target genes of the TF
    fig = dc.pl.volcano(
        data=df_results,
        x="log2FoldChange",
        y="padj",
        net=collectri,
        name=top_tf,
        top=30,
        thr_stat=snakemake.params["log2_fold_change_threshold"],
        thr_sign=snakemake.params["significance_threshold"],
        figsize=(10, 7.5),
        return_fig=True,
    )
    fig.suptitle(f"Target genes of {top_tf}", fontsize=12)
    volcano_file = tf_plots_dir / f"{top_tf}_volcano.svg"
    fig.savefig(volcano_file, dpi=500, bbox_inches="tight")
    plt.close(fig)

# perform pathway enrichment analysis
progeny = pd.read_csv(snakemake.input["progeny"], index_col=0)
pathway_acts, pathway_pvals = dc.mt.ulm(
    data=mat,
    net=progeny,
    tmin=snakemake.params["pathway_overlap_count"],
)

fig = dc.pl.barplot(
    data=pathway_acts,
    name=treated_vs_untreated_identifier,
    top=snakemake.params["top_pathways"],
    vertical=False,
    figsize=(7, 3),
    return_fig=True,
)
fig.tight_layout(pad=1)
pathways_barplot_path = Path(snakemake.output["pathways_barplot"])
fig.savefig(pathways_barplot_path, dpi=500, bbox_inches="tight")
plt.close(fig)
pathway_acts.to_csv(snakemake.output["pathways_table"])

# Create leading edge plots for top pathways
pathway_acts_ranked = pathway_acts.T.copy()
pathway_acts_ranked["activity_absolute"] = np.abs(
    pathway_acts_ranked[treated_vs_untreated_identifier]
)
top_pathways = (
    pathway_acts_ranked.sort_values("activity_absolute", ascending=False)
    .head(snakemake.params["top_pathways"])
    .index
)

pathways_plots_dir = pathways_barplot_path.parent

for top_pathway in top_pathways:
    # Get positive and negative edges if they exist
    pos_net = progeny[(progeny["source"] == top_pathway) & (progeny["weight"] > 0)]
    neg_net = progeny[(progeny["source"] == top_pathway) & (progeny["weight"] < 0)]

    # Create leading edge plots for positive and negative components
    if len(pos_net) > 0:
        try:
            _, pos_le = dc.pl.leading_edge(
                df_results,
                stat="stat",
                net=pos_net,
                name=top_pathway,
                figsize=(10, 5),
                return_fig=True,
            )
            plt.savefig(
                pathways_plots_dir / f"{top_pathway}_positive_leading_edge.svg",
                dpi=500,
                bbox_inches="tight",
            )
            plt.close()
        except Exception as e:
            warning(
                f"Could not create positive leading edge plot for {top_pathway}: {e}"
            )

    if len(neg_net) > 0:
        try:
            _, neg_le = dc.pl.leading_edge(
                df_results,
                stat="stat",
                net=neg_net,
                name=top_pathway,
                figsize=(10, 5),
                return_fig=True,
            )
            plt.savefig(
                pathways_plots_dir / f"{top_pathway}_negative_leading_edge.svg",
                dpi=500,
                bbox_inches="tight",
            )
            plt.close()
        except Exception as e:
            warning(
                f"Could not create negative leading edge plot for {top_pathway}: {e}"
            )

# Save processed df_results for downstream geneset analysis
df_results.to_csv(snakemake.output["processed_results"])
