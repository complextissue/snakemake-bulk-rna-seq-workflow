"""Create a transcript to gene map from a GTF annotation file."""

from pytximport.utils import create_transcript_to_gene_map
from snakemake.script import snakemake

df_transcript_gene_map = create_transcript_to_gene_map(
    species="human",
    host="may2024.archive.ensembl.org",
    target_field="gene_name",
)

df_transcript_gene_map.to_csv(
    snakemake.output[0],
    sep="\t",
    index=False,
)
