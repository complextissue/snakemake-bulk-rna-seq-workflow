# Helper functions for dynamic input based on read type
def get_salmon_input_reads(wildcards):
    """Get salmon input reads based on paired/unpaired detection."""
    if IS_PAIRED:
        return {
            "r1": f"resources/reads/trimmed/{wildcards.sample_id}_1.fastq.gz",
            "r2": f"resources/reads/trimmed/{wildcards.sample_id}_2.fastq.gz",
        }
    else:
        return {
            "r": f"resources/reads/trimmed/{wildcards.sample_id}.fastq.gz",
        }


# Salmon index files (shared between paired and unpaired)
SALMON_INDEX = multiext(
    "resources/reference/salmon/",
    "complete_ref_lens.bin",
    "ctable.bin",
    "ctg_offsets.bin",
    "duplicate_clusters.tsv",
    "info.json",
    "mphf.bin",
    "pos.bin",
    "pre_indexing.log",
    "rank.bin",
    "refAccumLengths.bin",
    "ref_indexing.log",
    "reflengths.bin",
    "refseq.bin",
    "seq.bin",
    "versionInfo.json",
)


rule quantify_reads_salmon:
    input:
        unpack(get_salmon_input_reads),
        index=SALMON_INDEX,
    output:
        quant="resources/reads/quantified_salmon/{sample_id}/quant.sf",
        lib="resources/reads/quantified_salmon/{sample_id}/lib_format_counts.json",
        info="resources/reads/quantified_salmon/{sample_id}/aux_info/meta_info.json",
        flenDist="resources/reads/quantified_salmon/{sample_id}/libParams/flenDist.txt",
    log:
        "results/logs/quantify_reads_salmon/{sample_id}.log",
    benchmark:
        "results/benchmarks/quantify_reads_salmon/{sample_id}.tsv"
    params:
        libtype=config["quantify_reads_salmon"]["libtype"],
        extra=config["quantify_reads_salmon"]["extra"],
    threads: config["quantify_reads_salmon"]["threads"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000,  # 8GB, scales with retries
        runtime=lambda wildcards, attempt: attempt * 60,  # 60 min, scales with retries
    wrapper:
        "v7.8.1/bio/salmon/quant"
