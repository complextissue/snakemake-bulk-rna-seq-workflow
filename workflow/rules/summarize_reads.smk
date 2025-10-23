rule run_pytximport_salmon:
    input:
        files=expand(
            "resources/reads/quantified_salmon/{sample_id}/quant.sf",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_pytximport/counts_salmon_{counts_from_abundance_pytximport}.{output_format}",
    params:
        data_type="salmon",
        gene_level=False,
        inferential_replicates=True,
        counts_from_abundance="{counts_from_abundance_pytximport}",
    threads: config["summarize_reads"]["threads"]
    conda:
        "../envs/summarize_reads.yaml"
    script:
        "../scripts/summarize_reads.py"
