rule get_transcriptome:
    output:
        "resources/reference/transcriptome.fasta",
    params:
        species=config["get_reference"]["species"],
        datatype="cdna",
        build=config["get_reference"]["build"],
        release=config["get_reference"]["release"],
    log:
        "results/logs/get_reference/get_transcriptome.log",
    cache: "omit-software"
    wrapper:
        "v3.13.6/bio/reference/ensembl-sequence"


rule get_genome:
    output:
        "resources/reference/genome.fasta",
    params:
        species=config["get_reference"]["species"],
        datatype="dna",
        build=config["get_reference"]["build"],
        release=config["get_reference"]["release"],
    log:
        "results/logs/get_reference/get_genome.log",
    cache: "omit-software"
    wrapper:
        "v3.13.6/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/reference/annotation.gtf",
    params:
        species=config["get_reference"]["species"],
        build=config["get_reference"]["build"],
        release=config["get_reference"]["release"],
    log:
        "results/logs/get_reference/get_annotation.log",
    cache: "omit-software"
    wrapper:
        "v3.13.6/bio/reference/ensembl-annotation"


rule create_transcript_to_gene_map:
    input:
        gtf="resources/reference/annotation.gtf",
    output:
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    conda:
        "../envs/get_reference.yaml"
    script:
        "../scripts/create_transcript_to_gene_map.py"
