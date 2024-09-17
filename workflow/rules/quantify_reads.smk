rule quantify_reads_salmon:
    input:
        r1="resources/reads/trimmed/{sample_id}_1.fastq.gz",
        r2="resources/reads/trimmed/{sample_id}_2.fastq.gz",
        index=multiext(
            "resources/reference/transcriptome_index/",
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


rule quantify_reads_bowtie:
    input:
        sample=["resources/reads/trimmed/{sample_id}_1.fastq.gz", "resources/reads/trimmed/{sample_id}_2.fastq.gz"],
        idx=multiext(
            "resources/reference/bowtie/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        "resources/reads/quantified_bowtie/{sample_id}.bam",
    log:
        "results/logs/quantify_reads_bowtie/{sample_id}.log",
    params:
        extra="",
    threads: config["quantify_reads_bowtie"]["threads"]
    wrapper:
        "v4.3.0/bio/bowtie2/align"


rule quantify_reads_rsem:
    input:
        bam="resources/reads/quantified_bowtie/{sample_id}.bam",
        reference=multiext(
            "resources/reference/rsem/reference", ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"
        ),
    output:
        genes_results="resources/reads/quantified_rsem/{sample_id}.genes.results",
        isoforms_results="resources/reads/quantified_rsem/{sample_id}.isoforms.results",
    log:
        "results/logs/quantify_reads_rsem/{sample_id}.log",
    params:
        paired_end=True,
        extra="--seed 42",
    threads: config["quantify_reads_rsem"]["threads"]
    wrapper:
        "v4.3.0/bio/rsem/calculate-expression"


rule quantify_reads_star:
    input:
        fq1="resources/reads/{accession}_1.fastq.gz",
        fq2="resources/reads/{accession}_2.fastq.gz",
        idx="resources/reference/star/index/",
    output:
        output_path="resources/reads/quantified_star/{accession}/",
        aln="resources/reads/quantified_star/{accession}/aligned.bam",
        aln_transcriptome="resources/reads/quantified_star/{accession}/Aligned.toTranscriptome.out.bam",
    log:
        "results/logs/quantify_reads_star/{accession}.log",
    params:
        extra=(
            "--outSAMtype BAM Unsorted "
            "--quantMode TranscriptomeSAM GeneCounts "
            "--outFilterScoreMin 30 "
            "--outFilterMultimapNmax 10 --winAnchorMultimapNmax 50 "
            "--alignEndsType Local "
        ),
    threads: config["quantify_reads_star"]["threads"]
    conda:
        "../envs/quantify_reads_star.yaml"
    script:
        "../scripts/quantify_reads_star.py"


rule quantify_reads_rsem_star:
    input:
        bam="resources/reads/quantified_star/{sample_id}/Aligned.toTranscriptome.out.bam",
        reference=multiext(
            "resources/reference/rsem/reference", ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"
        ),
    output:
        genes_results="resources/reads/quantified_rsem_star/{sample_id}.genes.results",
        isoforms_results="resources/reads/quantified_rsem_star/{sample_id}.isoforms.results",
    log:
        "results/logs/quantify_reads_rsem_star/{sample_id}.log",
    params:
        paired_end=True,
        extra="--seed 42",
    threads: config["quantify_reads_rsem"]["threads"]
    wrapper:
        "v4.3.0/bio/rsem/calculate-expression"


rule quantify_reads_kallisto:
    input:
        fastq=["resources/reads/trimmed/{sample_id}_1.fastq.gz", "resources/reads/trimmed/{sample_id}_2.fastq.gz"],
        index="resources/reference/kallisto/transcriptome.idx",
    output:
        directory("resources/reads/quantified_kallisto/{sample_id}"),
    params:
        extra="",
    log:
        "results/logs/quantify_reads_kallisto/{sample_id}.log",
    threads: config["quantify_reads_kallisto"]["threads"]
    wrapper:
        "v4.3.0/bio/kallisto/quant"
