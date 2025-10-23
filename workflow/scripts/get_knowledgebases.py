"""Perform gene set enrichment analysis on the results of a differential expression analysis using decoupler-py."""

import decoupler as dc
from snakemake.script import snakemake

collectri = dc.op.collectri(
    organism=snakemake.params["species"],
    remove_complexes=False,
)
collectri.to_csv(snakemake.output["collectri"])

progeny = dc.op.progeny(
    organism=snakemake.params["species"],
    top=500,
)
progeny.to_csv(snakemake.output["progeny"])

msigdb = dc.op.resource(
    "MSigDB",
    organism=snakemake.params["species"],
)
# Filter to keep only rows where genesymbol is present
msigdb = msigdb[~msigdb.duplicated(["geneset", "genesymbol"])]
msigdb = msigdb.dropna(subset=["genesymbol"])
msigdb.to_csv(snakemake.output["msigdb"])
