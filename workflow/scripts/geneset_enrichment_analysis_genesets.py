"""Perform gene set enrichment analysis for MSigDB genesets."""

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

# Prepare data matrix: transpose stat column so genes are columns
mat = df_results[["stat"]].T

# Perform GSEA using ULM (Univariate Linear Model)
enriched_genesets, msigdb_padj = dc.mt.ulm(
    data=mat,
    net=msigdb_renamed,
    tmin=3,
)

# Filter by adjusted p-value threshold (< 0.05)
msk = (msigdb_padj.T < 0.05).iloc[:, 0]
enriched_genesets = enriched_genesets.loc[:, msk]
msigdb_padj = msigdb_padj.loc[:, msk]

# Create barplot only if we have significant results
if len(enriched_genesets.columns) > 0:
    fig = dc.pl.barplot(
        data=enriched_genesets,
        name="stat",
        top=snakemake.params["top_genesets"],
        vertical=True,
        figsize=(15, 10),
        return_fig=True,
    )
    fig.tight_layout(pad=1)
else:
    # Create empty plot if no significant results
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.text(
        0.5,
        0.5,
        "No significantly enriched genesets found (padj < 0.05)",
        ha="center",
        va="center",
        fontsize=12,
        transform=ax.transAxes,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.tight_layout(pad=1)

geneset_barplot_path = Path(snakemake.output["geneset_dotplot"])
fig.savefig(geneset_barplot_path, dpi=500, bbox_inches="tight")
plt.close(fig)

# Save enrichment results with p-values
enrichment_results = enriched_genesets.T.copy()
enrichment_results["p_value"] = msigdb_padj.T.values
enrichment_results.to_csv(snakemake.output["geneset_table"])
