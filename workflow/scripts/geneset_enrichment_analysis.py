"""Perform gene set enrichment analysis on the results of a differential expression analysis using decoupler-py."""

from logging import warning
from pathlib import Path

import decoupler as dc
import numpy as np
import pandas as pd
from snakemake.script import snakemake

df_results = pd.read_csv(snakemake.input["results_table"], index_col=0, header=0)

# replace gene ids with gene names
if "gene_name" in df_results.columns:
    if df_results["gene_name"].nunique() != df_results.index.nunique():
        warning(
            "Missing or duplicated gene names, setting to gene_ids:\n"
            f"{df_results[df_results.duplicated('gene_name')].head(5)}"
        )
        df_results.loc[df_results.duplicated("gene_name"), "gene_name"] = df_results[
            df_results.duplicated("gene_name")
        ].index

    df_results = df_results.set_index("gene_name")

# save a volcano plot of the differentially expressed genes
fig = dc.plot_volcano_df(
    df_results,
    x="log2FoldChange",
    y="padj",
    sign_thr=snakemake.params["significance_threshold"],
    lFCs_thr=snakemake.params["log2_fold_change_threshold"],
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
fig.suptitle(f"Differentially expressed genes ({count_up} up, {count_down} down)", fontsize=12)
fig.savefig(snakemake.output["volcano_plot"], dpi=500)

# perform transcription factor enrichment analysis
collectri = dc.get_collectri(
    organism=snakemake.params["species"],
    split_complexes=False,
)

treated_vs_untreated_identifier = f"{snakemake.params['treated_name']}.vs.{snakemake.params['untreated_name']}"

mat = df_results[["stat"]].T.rename(
    index={"stat": treated_vs_untreated_identifier},
)

tf_acts, tf_pvals = dc.run_ulm(
    mat=mat,
    net=collectri,
    verbose=True,
    min_n=10,
)
fig = dc.plot_barplot(
    acts=tf_acts,
    contrast=treated_vs_untreated_identifier,
    top=snakemake.params["top_transcription_factors"],
    vertical=False,
    figsize=(10, 5),
    return_fig=True,
)
fig.tight_layout(pad=1)
transcription_factors_barplot_path = Path(snakemake.output["transcription_factors_barplot"])
fig.savefig(transcription_factors_barplot_path, dpi=500)
tf_acts.to_csv(snakemake.output["transcription_factors_table"])

# plot the log2 fold change of the target genes of the top transcription factors by their absolute log2 fold change
tf_acts = tf_acts.T
tf_acts["activity_absolute"] = np.abs(tf_acts[treated_vs_untreated_identifier])
top_tfs = (
    tf_acts.sort_values("activity_absolute", ascending=False).head(snakemake.params["top_transcription_factors"]).index
)

# Extract logFCs and pvals
log2_fold_changes = df_results[["log2FoldChange"]].T.rename(index={"log2FoldChange": treated_vs_untreated_identifier})
p_values_adjusted = df_results[["padj"]].T.rename(index={"padj": treated_vs_untreated_identifier})

for top_tf in top_tfs:
    fig = dc.plot_volcano(
        log2_fold_changes,
        p_values_adjusted,
        treated_vs_untreated_identifier,
        name=top_tf,
        net=collectri,
        top=30,
        sign_thr=snakemake.params["significance_threshold"],
        lFCs_thr=snakemake.params["log2_fold_change_threshold"],
        figsize=(10, 7.5),
        return_fig=True,
    )
    fig.suptitle(f"Target genes of {top_tf}", fontsize=12)
    fig.savefig(
        transcription_factors_barplot_path.parent / f"{top_tf}{transcription_factors_barplot_path.suffix}", dpi=500
    )

# perform pathway enrichment analysis
progeny = dc.get_progeny(
    organism=snakemake.params["species"],
    top=500,
)
pathway_acts, pathway_pvals = dc.run_mlm(
    mat=mat,
    net=progeny,
    min_n=snakemake.params["pathway_overlap_count"],
)
fig = dc.plot_barplot(
    acts=pathway_acts,
    contrast=treated_vs_untreated_identifier,
    top=snakemake.params["top_pathways"],
    vertical=False,
    figsize=(7, 3),
    return_fig=True,
)
fig.tight_layout(pad=1)
fig.savefig(snakemake.output["pathways_barplot"], dpi=500)
pathway_acts.to_csv(snakemake.output["pathways_table"])

# perform geneset enrichment analysis on the Reactome resource
msigdb_all = dc.get_resource(
    "MSigDB",
    organism=snakemake.params["species"],
)
msigdb = msigdb_all.copy()[msigdb_all["collection"] == snakemake.params["msigdb_geneset"]]
msigdb = msigdb[~msigdb.duplicated(["geneset", "genesymbol"])]
msigdb = msigdb.dropna()

if msigdb["geneset"].str.contains("REACTOME_").any():
    msigdb["geneset"] = msigdb["geneset"].apply(lambda geneset_name: geneset_name.replace("REACTOME_", ""))

if msigdb["geneset"].str.contains("KEGG_").any():
    msigdb["geneset"] = msigdb["geneset"].apply(lambda geneset_name: geneset_name.replace("KEGG_", ""))

if msigdb["geneset"].str.contains("GOBP_").any():
    msigdb["geneset"] = msigdb["geneset"].apply(lambda geneset_name: geneset_name.replace("GOBP_", ""))

msigdb["geneset"] = msigdb["geneset"].apply(lambda geneset_name: geneset_name.replace("_", " "))
top_genes = df_results[df_results["padj"] < snakemake.params["significance_threshold"]]

enriched_genesets = dc.get_ora_df(
    df=top_genes,
    net=msigdb,
    source="geneset",
    target="genesymbol",
)
fig = dc.plot_dotplot(
    enriched_genesets.sort_values("Combined score", ascending=False).head(snakemake.params["top_genesets"]),
    x="Combined score",
    y="Term",
    s="Odds ratio",
    c="FDR p-value",
    scale=0.75,
    figsize=(12, 5),
    return_fig=True,
)
fig.tight_layout(pad=1)
fig.savefig(snakemake.output["geneset_dotplot"], dpi=500)
enriched_genesets.to_csv(snakemake.output["geneset_table"])
