from snakemake.utils import validate
import pandas as pd


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

#
# Load samples
#

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample_id", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

reads = ["1", "2"]


wildcard_constraints:
    sample_id="|".join(samples.index),
    read="|".join(reads),


def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """
    wanted_input = []

    if config["check_quality"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/{sample_id}_{read}.fastq.gz",
                    "results/plots/check_quality/run_fastqc/{sample_id}_{read}.html",
                    "results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
                    "results/logs/check_quality/run_fastqc/{sample_id}_{read}.log",
                    "results/plots/check_quality/multiqc/report.html",
                    "results/plots/check_quality/run_fastp/{sample_id}.html",
                ],
                sample_id=samples.index,
                read=reads,
            )
        )

    if config["get_reference"]["run"]:
        wanted_input.extend(
            [
                "resources/reference/transcriptome.fasta",
                "resources/reference/genome.fasta",
                "resources/reference/annotation.gtf",
                "resources/reference/transcript_to_gene_map.tsv",
            ]
        )

    if config["build_index"]["run"]:
        wanted_input.extend(
            [
                "resources/reference/gentrome.fasta",
                "resources/reference/decoys.txt",
                "resources/reference/transcriptome_index/info.json",
                "results/logs/build_index/create_salmon_decoys.log",
                "results/logs/build_index/create_salmon_index.log",
            ]
        )

    if config["quantify_reads"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/quantified/{sample_id}/quant.sf",
                    "resources/reads/quantified/{sample_id}/lib_format_counts.json",
                    "results/logs/quantify_reads/{sample_id}.log",
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
            [
                "resources/reads/summarized/counts.h5ad",
            ],
        )

    if config["summarize_reads"]["run_tximport"]:
        wanted_input.extend(
            [
                "results/tables/summarize_reads_tximport/counts.csv",
                "results/tables/perform_dge_analysis/pydeseq2_tximport.csv",
                "results/plots/perform_dge_analysis/pca_plot_tximport.svg",
            ],
        )

    if config["perform_dge_analysis"]["run"]:
        wanted_input.extend(
            [
                "results/tables/perform_dge_analysis/pydeseq2.csv",
                "results/plots/perform_dge_analysis/pca_plot.svg",
            ],
        )

    if config["perform_gse_analysis"]["run"]:
        wanted_input.extend(
            [
                "results/plots/perform_dge_analysis/volcano_plot.svg",
                "results/plots/perform_gse_analysis/collectri_barplot.svg",
                "results/tables/perform_gse_analysis/collectri.csv",
                "results/plots/perform_gse_analysis/progeny_barplot.svg",
                "results/tables/perform_gse_analysis/progeny.csv",
                f"results/plots/perform_gse_analysis/{config['perform_gse_analysis']['msigdb_geneset']}_dotplot.svg",
                f"results/tables/perform_gse_analysis/{config['perform_gse_analysis']['msigdb_geneset']}.csv",
            ],
        )

    return wanted_input
