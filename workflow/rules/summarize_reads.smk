rule run_pytximport:
    input:
        files=expand(
            "resources/reads/quantified_salmon/{sample_id}/quant.sf",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        adata=f"resources/reads/summarized_pytximport/counts_{config['summarize_reads']['counts_from_abundance']}.h5ad",
        counts=f"resources/reads/summarized_pytximport/counts_{config['summarize_reads']['counts_from_abundance']}.csv",
    params:
        return_transcripts=config["summarize_reads"]["return_transcripts"],
        counts_from_abundance=config["summarize_reads"]["counts_from_abundance"],
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
        counts=f"resources/reads/summarized_tximport/counts_{config['summarize_reads_tximport']['counts_from_abundance']}.csv",
    params:
        counts_from_abundance=config["summarize_reads_tximport"]["counts_from_abundance"],
        return_transcripts=config["summarize_reads_tximport"]["return_transcripts"],
    conda:
        "../envs/summarize_reads_tximport.yaml"
    script:
        "../scripts/summarize_reads_tximport.R"
