rule get_knowledgebases:
    output:
        collectri="results/tables/perform_gse_analysis/collectri.csv",
        progeny="results/tables/perform_gse_analysis/progeny.csv",
        msigdb="results/tables/perform_gse_analysis/msigdb.csv",
    params:
        species=config["experiment"]["species"],
    conda:
        "../envs/perform_gse_analysis.yaml"
    script:
        "../scripts/get_knowledgebases.py"


rule run_decoupler_base:
    input:
        results_table="results/tables/perform_dge_analysis/pydeseq2.csv",
        collectri="results/tables/perform_gse_analysis/collectri.csv",
        progeny="results/tables/perform_gse_analysis/progeny.csv",
    output:
        volcano_plot="results/plots/perform_dge_analysis/volcano_plot_base.svg",
        transcription_factors_barplot="results/plots/perform_gse_analysis/collectri_barplot_base.svg",
        transcription_factors_table="results/tables/perform_gse_analysis/collectri_base.csv",
        pathways_barplot="results/plots/perform_gse_analysis/progeny_barplot_base.svg",
        pathways_table="results/tables/perform_gse_analysis/progeny_base.csv",
        processed_results="results/tables/perform_gse_analysis/processed_results.csv",
    params:
        treated_name=config["experiment"]["treated_name"],
        untreated_name=config["experiment"]["untreated_name"],
        significance_threshold=config["perform_gse_analysis"]["significance_threshold"],
        log2_fold_change_threshold=config["perform_gse_analysis"][
            "log2_fold_change_threshold"
        ],
        top_genes=config["perform_gse_analysis"]["top_genes"],
        top_transcription_factors=config["perform_gse_analysis"][
            "top_transcription_factors"
        ],
        pathway_overlap_count=config["perform_gse_analysis"]["pathway_overlap_count"],
        top_pathways=config["perform_gse_analysis"]["top_pathways"],
    conda:
        "../envs/perform_gse_analysis.yaml"
    script:
        "../scripts/geneset_enrichment_analysis_base.py"


rule run_decoupler_genesets:
    input:
        processed_results="results/tables/perform_gse_analysis/processed_results.csv",
        msigdb="results/tables/perform_gse_analysis/msigdb.csv",
    output:
        geneset_dotplot="results/plots/perform_gse_analysis/{msigdb_geneset}_dotplot.svg",
        geneset_table="results/tables/perform_gse_analysis/{msigdb_geneset}.csv",
    wildcard_constraints:
        msigdb_geneset="|".join(config["perform_gse_analysis"]["msigdb_geneset"]),
    params:
        species=config["experiment"]["species"],
        significance_threshold=config["perform_gse_analysis"]["significance_threshold"],
        msigdb_geneset="{msigdb_geneset}",
        top_genesets=config["perform_gse_analysis"]["top_genesets"],
    conda:
        "../envs/perform_gse_analysis.yaml"
    script:
        "../scripts/geneset_enrichment_analysis_genesets.py"
