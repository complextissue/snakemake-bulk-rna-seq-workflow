rule run_pytximport:
    input:
        files=expand(
            "resources/reads/quantified/{sample_id}/quant.sf",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized/counts.h5ad",
    params:
        return_transcripts=config["summarize_reads"]["return_transcripts"],
        counts_from_abundance=config["summarize_reads"]["counts_from_abundance"],
    conda:
        "../envs/summarize_reads.yaml"
    script:
        "../scripts/summarize_reads.py"
