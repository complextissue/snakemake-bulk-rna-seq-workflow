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

txi <- tximport(
    files,
    type = "salmon",
    txOut = tx_out,
    tx2gene = tx2gene,
    countsFromAbundance = snakemake@params[["counts_from_abundance"]],
    dropInfReps = FALSE,
    infRepStat = rowMedians,
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE
)

# export the counts to a .csv file
write.csv(txi$counts, file = snakemake@output[["counts"]])
