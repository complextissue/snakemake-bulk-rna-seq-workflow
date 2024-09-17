rule run_pytximport:
    input:
        files=expand(
            "resources/reads/quantified_salmon/{sample_id}/quant.sf",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_pytximport/counts_{counts_from_abundance_pytximport}.{output_format}",
    params:
        counts_from_abundance="{counts_from_abundance_pytximport}",
    conda:
        "../envs/summarize_reads.yaml"
    script:
        "../scripts/summarize_reads.py"


rule run_tximport:
    input:
        files=expand(
            "resources/reads/quantified_salmon/{sample_id}/quant.sf",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_tximport/counts_{counts_from_abundance_tximport}.csv",
    params:
        counts_from_abundance="{counts_from_abundance_tximport}",
    conda:
        "../envs/summarize_reads_tximport.yaml"
    script:
        "../scripts/summarize_reads_tximport.R"
