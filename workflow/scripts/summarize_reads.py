"""Summarize transcript-level read counts to gene-level counts using pytximport."""

from pathlib import Path

import numpy as np
from pytximport import tximport
from snakemake.script import snakemake

output_path = Path(snakemake.output[0])

if output_path.suffix == ".h5ad":
    output_format = "h5ad"
elif output_path.suffix == ".csv":
    output_format = "csv"
else:
    raise ValueError(f"Unknown file format: {output_path.suffix}. Supported formats are .h5ad and .csv.")

config = {
    "file_paths": snakemake.input["files"],
    "data_type": snakemake.params["data_type"],
    "gene_level": snakemake.params["gene_level"],
    "transcript_gene_map": snakemake.input["transcript_to_gene_map"],
    "output_path": snakemake.output[0],
    "return_transcript_data": snakemake.params["counts_from_abundance"] not in ["None", "length_scaled_tpm"],
    "counts_from_abundance": None if snakemake.params["counts_from_abundance"] == "None" else snakemake.params["counts_from_abundance"],
    "inferential_replicates": snakemake.params["inferential_replicates"],
    "inferential_replicate_variance":snakemake.params["inferential_replicates"],
    "inferential_replicate_transformer": lambda x: np.median(x, axis=1) if snakemake.params["inferential_replicates"] else None,
    "output_format": output_format,
}

result = tximport(
    file_paths=config["file_paths"],
    data_type=config["data_type"],
    transcript_gene_map=config["transcript_gene_map"],
    return_transcript_data=config["return_transcript_data"],
    gene_level=config["gene_level"],
    inferential_replicates=config["inferential_replicates"],
    inferential_replicate_variance=config["inferential_replicate_variance"],
    inferential_replicate_transformer=config["inferential_replicate_transformer"],
    counts_from_abundance=config["counts_from_abundance"],
    ignore_after_bar=True,
    ignore_transcript_version=True,
    output_type="anndata",
    output_format=config["output_format"],
    output_path=config["output_path"],
)
