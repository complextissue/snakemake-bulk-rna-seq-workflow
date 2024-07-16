library(tximport)
library(readr)

tx2gene <- read_tsv(
    snakemake@input[['transcript_to_gene_map']],
    show_col_types = FALSE,
)

rowMedians <- function(x) {
    apply(x, 1, median, na.rm = TRUE)
}

files <- sapply(
    snakemake@input[["files"]],
    function(file) as.character(file)
);

txi <- tximport(
    files,
    type = "salmon",
    txOut = snakemake@params[['return_transcripts']],
    tx2gene = tx2gene,
    countsFromAbundance = snakemake@params[['counts_from_abundance']],
    dropInfReps = FALSE,
    infRepStat = rowMedians,
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE
)

# export the counts to a .csv file
write.csv(txi$counts, file = snakemake@output[['counts']])
