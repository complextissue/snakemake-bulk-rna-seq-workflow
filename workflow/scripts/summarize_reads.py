"""Summarize transcript-level read counts to gene-level counts using pytximport."""

from pathlib import Path

import numpy as np
from pytximport import tximport
from snakemake.script import snakemake

save_path = Path(snakemake.output[0])

if save_path.suffix == ".h5ad":
    output_format = "h5ad"
elif save_path.suffix == ".csv":
    output_format = "csv"
else:
    raise ValueError(f"Unknown file format: {save_path.suffix}. Supported formats are .h5ad and .csv.")

config = {
    "file_paths": snakemake.input["files"],
    "transcript_gene_map": snakemake.input["transcript_to_gene_map"],
    "save_path": snakemake.output[0],
    "return_transcript_data": snakemake.params["return_transcripts"],
    "counts_from_abundance": snakemake.params["counts_from_abundance"],
    "output_format": output_format,
}

result = tximport(
    file_paths=config["file_paths"],
    data_type="salmon",
    transcript_gene_map=config["transcript_gene_map"],
    return_transcript_data=config["return_transcript_data"],
    inferential_replicates=True,
    inferential_replicate_variance=True,
    inferential_replicate_transformer=lambda x: np.median(x, axis=1),
    counts_from_abundance=config["counts_from_abundance"],
    ignore_after_bar=True,
    ignore_transcript_version=True,
    output_type="anndata",
    output_format=config["output_format"],
    save_path=config["save_path"],
)
