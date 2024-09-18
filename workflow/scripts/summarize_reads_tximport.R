library(tximport)
library(readr)
library(matrixStats)

tx2gene <- read_tsv(
    snakemake@input[["transcript_to_gene_map"]],
    show_col_types = FALSE,
)

files <- sapply(
    snakemake@input[["files"]],
    function(file) as.character(file)
)

tx_out <- snakemake@params[["counts_from_abundance"]] %in% c(
    "scaledTPM",
    "dtuScaledTPM"
)

dropInfReps <- !snakemake@params[["inferential_replicates"]]
if (dropInfReps) {
    infRepStat <- rowMedians
    varReduce <- TRUE
} else {
    infRepStat <- NULL
    varReduce <- FALSE
}

txi <- tximport(
    files,
    type = snakemake@params[["data_type"]],
    txOut = tx_out,
    txIn = !snakemake@params[["gene_level"]],
    tx2gene = tx2gene,
    countsFromAbundance = snakemake@params[["counts_from_abundance"]],
    dropInfReps = dropInfReps,
    infRepStat = infRepStat,
    varReduce = varReduce,
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE
)

# export the counts to a .csv file
write.csv(txi$counts, file = snakemake@output[["counts"]])
