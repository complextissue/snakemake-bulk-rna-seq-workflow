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
        "v7.8.1/bio/salmon/decoys"


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
        "v7.8.1/bio/salmon/index"
