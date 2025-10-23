from snakemake.utils import validate
import pandas as pd
from pathlib import Path
import os


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

#
# Load samples
#

# Read metadata.csv and transform to match expected format
samples = pd.read_csv(config["samples"], sep=",")
samples = samples.set_index("sample_id", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

# Get paired vs single-end setting from config
IS_PAIRED = config["experiment"]["paired"]
reads = ["1", "2"] if IS_PAIRED else []

output_formats = ["h5ad", "csv"]
counts_from_abundances_pytximport = ["length_scaled_tpm"]

# Build wildcard constraints based on read type
if IS_PAIRED:

    wildcard_constraints:
        sample_id="|".join(samples.index),
        output_format="|".join(output_formats),
        counts_from_abundance_pytximport="|".join(counts_from_abundances_pytximport),
        read="|".join(reads),

else:

    wildcard_constraints:
        sample_id="|".join(samples.index),
        output_format="|".join(output_formats),
        counts_from_abundance_pytximport="|".join(counts_from_abundances_pytximport),


def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """
    wanted_input = []

    if config["check_quality"]["run"]:
        # Build base file list depending on read type
        if IS_PAIRED:
            base_files = [
                "resources/reads/raw/{sample_id}_{read}.fastq.gz",
                "results/plots/check_quality/run_fastp/{sample_id}.html",
            ]
            # Add FastQC files if enabled
            if config["check_quality"]["run_fastqc"]:
                base_files.extend(
                    [
                        "results/plots/check_quality/run_fastqc/{sample_id}_{read}.html",
                        "results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
                    ]
                )

            wanted_input.extend(
                expand(
                    base_files,
                    sample_id=samples.index,
                    read=reads,
                )
            )
        else:
            base_files = [
                "resources/reads/raw/{sample_id}.fastq.gz",
                "results/plots/check_quality/run_fastp/{sample_id}.html",
            ]
            # Add FastQC files if enabled
            if config["check_quality"]["run_fastqc"]:
                base_files.extend(
                    [
                        "results/plots/check_quality/run_fastqc/{sample_id}.html",
                        "results/plots/check_quality/run_fastqc/{sample_id}_fastqc.zip",
                    ]
                )

            wanted_input.extend(
                expand(
                    base_files,
                    sample_id=samples.index,
                )
            )

        # Add MultiQC report
        wanted_input.append("results/plots/check_quality/run_multiqc/report.html")

    if config["get_reference"]["run"]:
        wanted_input.extend(
            [
                "resources/reference/transcriptome.fasta",
                "resources/reference/genome.fasta",
                "resources/reference/annotation.gtf",
                "resources/reference/transcript_to_gene_map.tsv",
            ]
        )

    if config["build_index_salmon"]["run"]:
        wanted_input.extend(
            [
                "resources/reference/gentrome.fasta",
                "resources/reference/decoys.txt",
                "resources/reference/salmon/info.json",
                "results/logs/build_index/create_decoys_salmon.log",
                "results/logs/build_index/create_index_salmon.log",
            ]
        )

    if config["quantify_reads_salmon"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/quantified_salmon/{sample_id}/quant.sf",
                    "resources/reads/quantified_salmon/{sample_id}/lib_format_counts.json",
                    "results/logs/quantify_reads_salmon/{sample_id}.log",
                ],
                sample_id=samples.index,
            )
        )

        if config["check_quality"]["run"]:
            wanted_input.extend(
                [
                    "results/plots/check_quality/run_multiqc_after_salmon/report.html",
                ],
            )

    if config["summarize_reads"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/summarized_pytximport/counts_salmon_{counts_from_abundance_pytximport}.{output_format}",
                ],
                counts_from_abundance_pytximport=counts_from_abundances_pytximport,
                output_format=output_formats,
            )
        )

    if config["perform_dge_analysis"]["run"]:
        wanted_input.extend(
            [
                "results/tables/perform_dge_analysis/pydeseq2.csv",
                "results/plots/perform_dge_analysis/pca_plot.svg",
                "results/plots/perform_dge_analysis/pvalue_histogram.svg",
                "results/plots/perform_dge_analysis/ma_plot.svg",
            ],
        )

    if config["perform_gse_analysis"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "results/plots/perform_dge_analysis/volcano_plot_{msigdb_geneset}.svg",
                    "results/plots/perform_gse_analysis/collectri_barplot_{msigdb_geneset}.svg",
                    "results/tables/perform_gse_analysis/collectri_{msigdb_geneset}.csv",
                    "results/plots/perform_gse_analysis/progeny_barplot_{msigdb_geneset}.svg",
                    "results/tables/perform_gse_analysis/progeny_{msigdb_geneset}.csv",
                    "results/plots/perform_gse_analysis/{msigdb_geneset}_dotplot.svg",
                    "results/tables/perform_gse_analysis/{msigdb_geneset}.csv",
                ],
                msigdb_geneset=config["perform_gse_analysis"]["msigdb_geneset"],
            )
        )

    return wanted_input
