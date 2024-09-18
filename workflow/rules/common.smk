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
output_formats = ["h5ad", "csv"]
counts_from_abundances_pytximport = ["None", "length_scaled_tpm", "scaled_tpm", "dtu_scaled_tpm"]
counts_from_abundances_tximport = ["no", "lengthScaledTPM", "scaledTPM", "dtuScaledTPM"]


wildcard_constraints:
    sample_id="|".join(samples.index),
    read="|".join(reads),
    output_format="|".join(output_formats),
    counts_from_abundance_pytximport="|".join(counts_from_abundances_pytximport),
    counts_from_abundance_tximport="|".join(counts_from_abundances_tximport),


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
                    # "results/plots/check_quality/run_fastqc/{sample_id}_{read}.html",
                    # "results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
                    # "results/logs/check_quality/run_fastqc/{sample_id}_{read}.log",
                    "results/plots/check_quality/run_multiqc/report.html",
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

    if config["build_index_kallisto"]["run"]:
        wanted_input.extend(
            [
                "resources/reference/kallisto/transcriptome.idx",
                "results/logs/build_index/create_index_kallisto.log",
            ]
        )

    if config["build_index_rsem"]["run"]:
        wanted_input.extend(
            [
                directory("resources/reference/rsem/"),
                "resources/reference/rsem/reference.seq",
                *multiext(
                    "resources/reference/rsem/reference",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
                "results/logs/build_index/create_index_rsem.log",
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

    if config["quantify_reads_kallisto"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/quantified_kallisto/{sample_id}/",
                    # "resources/reads/quantified_kallisto/{sample_id}/abundance.h5",
                    # "resources/reads/quantified_kallisto/{sample_id}/run_info.json",
                    "results/logs/quantify_reads_kallisto/{sample_id}.log",
                ],
                sample_id=samples.index,
            )
        )

    if config["quantify_reads_rsem"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/quantified_rsem/{sample_id}.genes.results",
                    "resources/reads/quantified_rsem/{sample_id}.isoforms.results",
                    "results/logs/quantify_reads_rsem/{sample_id}.log",
                ],
                sample_id=samples.index,
            )
        )

    if config["summarize_reads"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/summarized_pytximport/counts_salmon_{counts_from_abundance_pytximport}.{output_format}",
                    "resources/reads/summarized_pytximport/counts_rsem_transcript_{counts_from_abundance_pytximport}.{output_format}",
                    "resources/reads/summarized_pytximport/counts_rsem_gene.{output_format}",
                    "resources/reads/summarized_pytximport/counts_kallisto_{counts_from_abundance_pytximport}.{output_format}",
                ],
                counts_from_abundance_pytximport=counts_from_abundances_pytximport,
                output_format=output_formats,
            )
        )

    if config["summarize_reads_tximport"]["run"]:
        wanted_input.extend(
            expand(
                [
                    "resources/reads/summarized_tximport/counts_salmon_{counts_from_abundance_tximport}.csv",
                    "resources/reads/summarized_tximport/counts_kallisto_{counts_from_abundance_tximport}.csv",
                    "resources/reads/summarized_tximport/counts_rsem_transcript_{counts_from_abundance_tximport}.csv",
                    "resources/reads/summarized_tximport/counts_rsem_gene.csv",
                ],
                counts_from_abundance_tximport=counts_from_abundances_tximport,
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
