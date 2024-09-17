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
        seq="resources/reference/rsem/reference.seq",
        grp="resources/reference/rsem/reference.grp",
        ti="resources/reference/rsem/reference.ti",
        other_1="resources/reference/rsem/reference.transcripts.fa",
        other_2="resources/reference/rsem/reference.idx.fa",
        other_3="resources/reference/rsem/reference.n2g.idx.fa",
    params:
        extra="--gtf resources/reference/rsem/annotations.gtf",
    log:
        "results/logs/build_index/create_index_rsem.log",
    threads: config["build_index_rsem"]["threads"]
    wrapper:
        "v4.3.0/bio/rsem/prepare-reference"


rule create_index_bowtie:
    input:
        ref="resources/reference/genome.fasta",
    output:
        multiext(
            "resources/reference/bowtie/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        extra="",
    log:
        "results/logs/build_index/create_index_bowtie.log",
    threads: config["build_index_bowtie"]["threads"]
    wrapper:
        "v4.3.0/bio/bowtie2/build"


rule create_index_star:
    input:
        fasta="resources/reference/genome.fasta",
        gtf="resources/reference/annotation.gtf",
    output:
        directory("resources/reference/star/index/"),
    message:
        "Testing STAR index"
    threads: config["build_index_star"]["threads"]
    params:
        sjdb_overhang=config["build_index_star"]["sjdb_overhang"],
        extra="",
    log:
        "results/logs/build_index_star/build_star_index.log",
    wrapper:
        "v4.3.0/bio/star/index"


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
