rule run_pydeseq2:
    input:
        counts="resources/reads/summarized/counts.h5ad",
        sample_table="config/samples.tsv",
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        results_table="results/tables/perform_dge_analysis/pydeseq2.csv",
        pca_plot="results/plots/perform_dge_analysis/pca_plot.svg",
    params:
        treated_name=config["experiment"]["treated_name"],
        untreated_name=config["experiment"]["untreated_name"],
        min_mean_counts=config["perform_dge_analysis"]["min_mean_counts"],
        min_max_counts=config["perform_dge_analysis"]["min_max_counts"],
        min_sample_counts=config["perform_dge_analysis"]["min_sample_counts"],
        refit_cooks=config["perform_dge_analysis"]["refit_cooks"],
        shrink_log2_fold_change=config["perform_dge_analysis"]["shrink_log2_fold_change"],
    threads: config["perform_dge_analysis"]["threads"]
    conda:
        "../envs/perform_dge_analysis.yaml"
    script:
        "../scripts/differential_gene_expression.py"


rule run_pydeseq2_tximport:
    input:
        counts="results/tables/summarize_reads_tximport/counts.csv",
        sample_table="config/samples.tsv",
        transcript_to_gene_map="resources/reference/transcript_to_gene_map.tsv",
    output:
        results_table="results/tables/perform_dge_analysis/pydeseq2_tximport.csv",
        pca_plot="results/plots/perform_dge_analysis/pca_plot_tximport.svg",
    params:
        min_mean_counts=config["perform_dge_analysis"]["min_mean_counts"],
        min_max_counts=config["perform_dge_analysis"]["min_max_counts"],
        min_sample_counts=config["perform_dge_analysis"]["min_sample_counts"],
        treated_name=config["experiment"]["treated_name"],
        untreated_name=config["experiment"]["untreated_name"],
        refit_cooks=config["perform_dge_analysis"]["refit_cooks"],
        shrink_log2_fold_change=config["perform_dge_analysis"]["shrink_log2_fold_change"],
    threads: config["perform_dge_analysis"]["threads"]
    conda:
        "../envs/perform_dge_analysis.yaml"
    script:
        "../scripts/differential_gene_expression.py"
