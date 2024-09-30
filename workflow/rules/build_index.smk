rule create_decoys_salmon:
    input:
        transcriptome="resources/reference/transcriptome.fasta",
        genome="resources/reference/genome.fasta",
    output:
        gentrome="resources/reference/gentrome.fasta",
        decoys="resources/reference/decoys.txt",
    threads: config["build_index_salmon"]["threads"]
    log:
        "results/logs/build_index/create_decoys_salmon.log",
    wrapper:
        "v4.3.0/bio/salmon/decoys"


rule create_index_salmon:
    input:
        sequences="resources/reference/gentrome.fasta",
        decoys="resources/reference/decoys.txt",
    output:
        multiext(
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
    params:
        extra="",
    log:
        "results/logs/build_index/create_index_salmon.log",
    threads: config["build_index_salmon"]["threads"]
    wrapper:
        "v4.3.0/bio/salmon/index"


rule create_index_rsem:
    input:
        reference_genome="resources/reference/genome.fasta",
    output:
        reference_name=directory("resources/reference/rsem/"),
        seq="resources/reference/rsem/reference.seq",
        grp="resources/reference/rsem/reference.grp",
        ti="resources/reference/rsem/reference.ti",
        other_1="resources/reference/rsem/reference.transcripts.fa",
        other_2="resources/reference/rsem/reference.idx.fa",
        other_3="resources/reference/rsem/reference.n2g.idx.fa",
        other_4="resources/reference/rsem/reference.chrlist",
        bowtie2_reference=multiext(
            "resources/reference/rsem/reference",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        extra="--gtf resources/reference/annotation.gtf",
    log:
        "results/logs/build_index/create_index_rsem.log",
    conda:
        "../envs/rsem.yaml"
    threads: config["build_index_rsem"]["threads"]
    shell:
        """
        rsem-prepare-reference --bowtie2 --num-threads {threads} {params.extra} {input.reference_genome} {output.reference_name}/reference > {log}
        """


rule create_index_kallisto:
    input:
        fasta="resources/reference/transcriptome.fasta",
    output:
        index="resources/reference/kallisto/transcriptome.idx",
    params:
        extra="",
    log:
        "results/logs/build_index/create_index_kallisto.log",
    threads: config["build_index_kallisto"]["threads"]
    wrapper:
        "v4.3.0/bio/kallisto/index"
