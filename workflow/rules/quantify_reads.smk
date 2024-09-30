rule quantify_reads_salmon:
    input:
        r1="resources/reads/trimmed/{sample_id}_1.fastq.gz",
        r2="resources/reads/trimmed/{sample_id}_2.fastq.gz",
        index=multiext(
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
        ),
    output:
        quant="resources/reads/quantified_salmon/{sample_id}/quant.sf",
        lib="resources/reads/quantified_salmon/{sample_id}/lib_format_counts.json",
        info="resources/reads/quantified_salmon/{sample_id}/aux_info/meta_info.json",
        flenDist="resources/reads/quantified_salmon/{sample_id}/libParams/flenDist.txt",
    log:
        "results/logs/quantify_reads_salmon/{sample_id}.log",
    params:
        libtype=config["quantify_reads_salmon"]["libtype"],
        extra=config["quantify_reads_salmon"]["extra"],
    threads: config["quantify_reads_salmon"]["threads"]
    wrapper:
        "v4.3.0/bio/salmon/quant"


rule quantify_reads_rsem:
    input:
        fq1="resources/reads/trimmed/{sample_id}_1.fastq.gz",
        fq2="resources/reads/trimmed/{sample_id}_2.fastq.gz",
        bowtie2_reference=multiext(
            "resources/reference/rsem/reference",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        genes_results="resources/reads/quantified_rsem/{sample_id}.genes.results",
        isoforms_results="resources/reads/quantified_rsem/{sample_id}.isoforms.results",
    log:
        "results/logs/quantify_reads_rsem/{sample_id}.log",
    params:
        output_prefix="resources/reads/quantified_rsem/{sample_id}",
        reference="resources/reference/rsem/reference",
        extra="--seed 42 --estimate-rspd",
    threads: config["quantify_reads_rsem"]["threads"]
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-calculate-expression --paired-end --bowtie2 --num-threads {threads} {params.extra} {input.fq1} {input.fq2} {params.reference} {params.output_prefix} > {log}
        """


rule quantify_reads_kallisto:
    input:
        fastq=["resources/reads/trimmed/{sample_id}_1.fastq.gz", "resources/reads/trimmed/{sample_id}_2.fastq.gz"],
        index="resources/reference/kallisto/transcriptome.idx",
    output:
        directory("resources/reads/quantified_kallisto/{sample_id}"),
    params:
        extra=config["quantify_reads_kallisto"]["extra"],
    log:
        "results/logs/quantify_reads_kallisto/{sample_id}.log",
    threads: config["quantify_reads_kallisto"]["threads"]
    wrapper:
        "v4.3.0/bio/kallisto/quant"
