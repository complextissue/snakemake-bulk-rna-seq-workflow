"""Perform gene set enrichment analysis on the results of a differential expression analysis using decoupler-py."""

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

treated_vs_untreated_identifier = (
    f"{snakemake.params['treated_name']}.vs.{snakemake.params['untreated_name']}"
)

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
pathway_acts, pathway_pvals = dc.mt.mlm(
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

# perform geneset enrichment analysis on the MSigDB resource
msigdb_all = pd.read_csv(snakemake.input["msigdb"], index_col=0)
msigdb = msigdb_all.copy()[
    msigdb_all["collection"] == snakemake.params["msigdb_geneset"]
]
msigdb = msigdb[~msigdb.duplicated(["geneset", "genesymbol"])]
msigdb = msigdb.dropna()

# Clean up geneset names by removing prefixes and replacing underscores
if msigdb["geneset"].str.contains("REACTOME_").any():
    msigdb["geneset"] = msigdb["geneset"].apply(
        lambda geneset_name: geneset_name.replace("REACTOME_", "")
    )

if msigdb["geneset"].str.contains("KEGG_").any():
    msigdb["geneset"] = msigdb["geneset"].apply(
        lambda geneset_name: geneset_name.replace("KEGG_", "")
    )

if msigdb["geneset"].str.contains("GOBP_").any():
    msigdb["geneset"] = msigdb["geneset"].apply(
        lambda geneset_name: geneset_name.replace("GOBP_", "")
    )

msigdb["geneset"] = msigdb["geneset"].apply(
    lambda geneset_name: geneset_name.replace("_", " ")
)

# Rename columns for decoupler compatibility
msigdb_renamed = msigdb.rename(columns={"geneset": "source", "genesymbol": "target"})

# For ORA, we only use significant genes
top_genes = df_results[df_results["padj"] < snakemake.params["significance_threshold"]]

# Diagnostic logging
print(f"Number of significant genes: {len(top_genes)}")
print(f"Sample significant gene names: {list(top_genes.index[:20])}")
print(f"Sample gene names from df_results: {list(df_results.index[:10])}")
print(f"Number of unique genes in msigdb: {msigdb_renamed['target'].nunique()}")
print(f"Sample gene names from msigdb: {list(msigdb_renamed['target'].unique()[:10])}")

# Check overlap with ALL df_results genes
shared_genes_all = set(msigdb_renamed["target"].unique()) & set(df_results.index)
print(
    f"Number of shared genes between msigdb and df_results (all genes): {len(shared_genes_all)}"
)

# Check overlap with SIGNIFICANT genes only
shared_genes_sig = set(msigdb_renamed["target"].unique()) & set(top_genes.index)
print(
    f"Number of shared genes between msigdb and top_genes (significant): {len(shared_genes_sig)}"
)
print(
    f"Sample shared significant genes: {list(shared_genes_sig)[:10] if shared_genes_sig else 'None'}"
)

# Try to run ORA, but handle cases where there's no overlap between gene names
try:
    # ORA expects a binary format where we pass significant genes
    enriched_genesets = dc.mt.ora(
        data=top_genes,
        net=msigdb_renamed,
        source="source",
        target="target",
        tmin=0,
    )
except (AssertionError, ValueError) as e:
    # If no shared genes, create an empty results dataframe
    warning(f"Could not perform ORA: {e}")
    warning(
        "This likely means not enough significant genes overlap with this MSigDB collection"
    )
    warning(f"MSigDB organism: {snakemake.params['species']}")
    warning(f"MSigDB collection: {snakemake.params['msigdb_geneset']}")
    warning(
        "Try: 1) Lower significance threshold, 2) Check gene name format, 3) Use different collection"
    )
    enriched_genesets = pd.DataFrame(
        {
            "Term": [],
            "Genes": [],
            "Combined score": [],
            "Odds ratio": [],
            "FDR p-value": [],
        }
    )

# Sort by combined score and get top genesets
enriched_genesets_sorted = (
    enriched_genesets.sort_values("Combined score", ascending=False)
    if len(enriched_genesets) > 0
    else enriched_genesets
)
top_genesets_data = enriched_genesets_sorted.head(snakemake.params["top_genesets"])

# Create dotplot only if we have results
if len(top_genesets_data) > 0:
    fig = dc.pl.dotplot(
        top_genesets_data,
        x="Combined score",
        y="Term",
        s="Odds ratio",
        c="FDR p-value",
        scale=0.75,
        figsize=(12, 5),
        return_fig=True,
    )
else:
    # Create empty plot if no results
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.text(
        0.5,
        0.5,
        "No enriched genesets found",
        ha="center",
        va="center",
        fontsize=12,
        transform=ax.transAxes,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

fig.tight_layout(pad=1)
geneset_dotplot_path = Path(snakemake.output["geneset_dotplot"])
fig.savefig(geneset_dotplot_path, dpi=500, bbox_inches="tight")
plt.close(fig)
enriched_genesets.to_csv(snakemake.output["geneset_table"])

# Create leading edge plots for top genesets
genesets_plots_dir = geneset_dotplot_path.parent
top_genesets = (
    enriched_genesets_sorted.head(snakemake.params["top_genesets"])["Term"].tolist()
    if len(enriched_genesets_sorted) > 0
    else []
)

for top_geneset in top_genesets:
    try:
        _, leading_edge_genes = dc.pl.leading_edge(
            df_results,
            stat="stat",
            net=msigdb_renamed,
            name=top_geneset,
            figsize=(10, 5),
            return_fig=True,
        )
        plt.savefig(
            genesets_plots_dir / f"{top_geneset}_leading_edge.svg",
            dpi=500,
            bbox_inches="tight",
        )
        plt.close()
    except Exception as e:
        warning(f"Could not create leading edge plot for {top_geneset}: {e}")
