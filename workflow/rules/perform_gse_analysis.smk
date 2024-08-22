rule run_decoupler:
    input:
        results_table="results/tables/perform_dge_analysis/pydeseq2.csv",
    output:
        volcano_plot="results/plots/perform_dge_analysis/volcano_plot.svg",
        transcription_factors_barplot="results/plots/perform_gse_analysis/collectri_barplot.svg",
        transcription_factors_table="results/tables/perform_gse_analysis/collectri.csv",
        pathways_barplot="results/plots/perform_gse_analysis/progeny_barplot.svg",
        pathways_table="results/tables/perform_gse_analysis/progeny.csv",
        geneset_dotplot=f"results/plots/perform_gse_analysis/{config['perform_gse_analysis']['msigdb_geneset']}_dotplot.svg",
        geneset_table=f"results/tables/perform_gse_analysis/{config['perform_gse_analysis']['msigdb_geneset']}.csv",
    params:
        species=config["experiment"]["species"],
        treated_name=config["experiment"]["treated_name"],
        untreated_name=config["experiment"]["untreated_name"],
        significance_threshold=config["perform_gse_analysis"]["significance_threshold"],
        log2_fold_change_threshold=config["perform_gse_analysis"]["log2_fold_change_threshold"],
        top_genes=config["perform_gse_analysis"]["top_genes"],
        top_transcription_factors=config["perform_gse_analysis"]["top_transcription_factors"],
        pathway_overlap_count=config["perform_gse_analysis"]["pathway_overlap_count"],
        top_pathways=config["perform_gse_analysis"]["top_pathways"],
        msigdb_geneset=config["perform_gse_analysis"]["msigdb_geneset"],
        top_genesets=config["perform_gse_analysis"]["top_genesets"],
    conda:
        "../envs/perform_gse_analysis.yaml"
    script:
        "../scripts/geneset_enrichment_analysis.py"
