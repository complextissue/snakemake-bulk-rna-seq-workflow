rule quantify_reads:
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
        quant="resources/reads/quantified/{sample_id}/quant.sf",
        lib="resources/reads/quantified/{sample_id}/lib_format_counts.json",
    log:
        "results/logs/quantify_reads/{sample_id}.log",
    params:
        libtype=config["quantify_reads"]["libtype"],
        extra=config["quantify_reads"]["extra"],
    threads: config["quantify_reads"]["threads"]
    wrapper:
        "v3.13.6/bio/salmon/quant"
