rule run_tximport:
    input:
        files=expand(
            "resources/reads/quantified/{sample_id}/quant.sf",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="results/tables/summarize_reads_tximport/counts.csv",
    params:
        counts_from_abundance=config["summarize_reads"]["counts_from_abundance_tximport"],
        return_transcripts=config["summarize_reads"]["return_transcripts"],
    conda:
        "../envs/summarize_reads_tximport.yaml"
    script:
        "../scripts/summarize_reads_tximport.R"
