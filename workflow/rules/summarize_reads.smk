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
    conda:
        "../envs/summarize_reads.yaml"
    script:
        "../scripts/summarize_reads.py"


rule run_pytximport_rsem_transcript:
    input:
        files=expand(
            "resources/reads/quantified_rsem/{sample_id}.isoforms.results",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_pytximport/counts_rsem_transcript_{counts_from_abundance_pytximport}.{output_format}",
    params:
        data_type="rsem",
        gene_level=False,
        inferential_replicates=False,
        counts_from_abundance="{counts_from_abundance_pytximport}",
    conda:
        "../envs/summarize_reads.yaml"
    script:
        "../scripts/summarize_reads.py"


rule run_pytximport_rsem_gene:
    input:
        files=expand(
            "resources/reads/quantified_rsem/{sample_id}.genes.results",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_pytximport/counts_rsem_gene.{output_format}",
    params:
        data_type="rsem",
        gene_level=False,
        inferential_replicates=False,
        counts_from_abundance="None",
    conda:
        "../envs/summarize_reads.yaml"
    script:
        "../scripts/summarize_reads.py"


rule run_pytximport_kallisto:
    input:
        files=expand(
            "resources/reads/quantified_kallisto/{sample_id}/abundance.h5",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_pytximport/counts_kallisto_{counts_from_abundance_pytximport}.{output_format}",
    params:
        data_type="kallisto",
        gene_level=False,
        inferential_replicates=False,
        counts_from_abundance="{counts_from_abundance_pytximport}",
    conda:
        "../envs/summarize_reads.yaml"
    script:
        "../scripts/summarize_reads.py"


rule run_tximport_salmon:
    input:
        files=expand(
            "resources/reads/quantified_salmon/{sample_id}/quant.sf",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_tximport/counts_salmon_{counts_from_abundance_tximport}.csv",
    params:
        data_type="salmon",
        gene_level=False,
        inferential_replicates=False,
        counts_from_abundance="{counts_from_abundance_tximport}",
    conda:
        "../envs/summarize_reads_tximport.yaml"
    script:
        "../scripts/summarize_reads_tximport.R"


rule run_tximport_kallisto:
    input:
        files=expand(
            "resources/reads/quantified_kallisto/{sample_id}/abundance.h5",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_tximport/counts_kallisto_{counts_from_abundance_tximport}.csv",
    params:
        data_type="kallisto",
        gene_level=False,
        inferential_replicates=False,
        counts_from_abundance="{counts_from_abundance_tximport}",
    conda:
        "../envs/summarize_reads_tximport.yaml"
    script:
        "../scripts/summarize_reads_tximport.R"


rule run_tximport_rsem_transcript:
    input:
        files=expand(
            "resources/reads/quantified_rsem/{sample_id}.isoforms.results",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_tximport/counts_rsem_transcript_{counts_from_abundance_tximport}.csv",
    params:
        data_type="rsem",
        gene_level=False,
        inferential_replicates=False,
        counts_from_abundance="{counts_from_abundance_tximport}",
    conda:
        "../envs/summarize_reads_tximport.yaml"
    script:
        "../scripts/summarize_reads_tximport.R"


rule run_tximport_rsem_gene:
    input:
        files=expand(
            "resources/reads/quantified_rsem/{sample_id}.genes.results",
            sample_id=samples.index,
        ),
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        counts="resources/reads/summarized_tximport/counts_rsem_gene.csv",
    params:
        data_type="rsem",
        gene_level=True,
        inferential_replicates=False,
        counts_from_abundance="no",
    conda:
        "../envs/summarize_reads_tximport.yaml"
    script:
        "../scripts/summarize_reads_tximport.R"
