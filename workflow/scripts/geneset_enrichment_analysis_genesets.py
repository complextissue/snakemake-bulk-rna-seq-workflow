"""Perform gene set enrichment analysis for MSigDB genesets."""

from logging import warning
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import pandas as pd
from snakemake.script import snakemake

# Load the processed results from the base analysis
df_results = pd.read_csv(snakemake.input["processed_results"], index_col=0, header=0)

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
