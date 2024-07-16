rule create_salmon_decoys:
    input:
        transcriptome="resources/reference/transcriptome.fasta",
        genome="resources/reference/genome.fasta",
    output:
        gentrome="resources/reference/gentrome.fasta",
        decoys="resources/reference/decoys.txt",
    threads: config["build_index"]["threads"]
    log:
        "results/logs/build_index/create_salmon_decoys.log",
    wrapper:
        "v3.13.6/bio/salmon/decoys"


rule create_salmon_index:
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
    log:
        "results/logs/build_index/create_salmon_index.log",
    threads: config["build_index"]["threads"]
    params:
        extra="",
    wrapper:
        "v3.13.6/bio/salmon/index"
