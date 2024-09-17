"""Create a transcript to gene map from a GTF annotation file."""

import pandas as pd
from pytximport.utils import create_transcript_to_gene_map
from snakemake.script import snakemake

df_transcript_gene_map = create_transcript_to_gene_map(
    species=snakemake.params["species"],
    host=snakemake.params["host"],
    target_field="external_gene_name",
)

df_transcript_gene_map_gene_id = create_transcript_to_gene_map(
    species=snakemake.params["species"],
    host=snakemake.params["host"],
    target_field="ensembl_gene_id",
)

# For the transcript_id values that are not in the first map, add them from the second map
df_transcript_gene_map_gene_id_only = df_transcript_gene_map_gene_id[
    ~df_transcript_gene_map_gene_id["transcript_id"].isin(
        df_transcript_gene_map["transcript_id"]
    )
]

df_transcript_gene_map = pd.concat([
    df_transcript_gene_map,
    df_transcript_gene_map_gene_id_only,
])

df_transcript_gene_map.to_csv(
    snakemake.output[0],
    sep="\t",
    index=False,
)
