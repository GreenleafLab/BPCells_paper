suppressPackageStartupMessages({
    library(BPCells)
    library(tibble)
    library(readr)
})


args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)
input_dir <- args[1]
output_dir <- args[2]
threads <- as.integer(args[3])

output_timing <- file.path(output_dir, "timing.tsv")

mat <- readRDS(file.path(input_dir, "matrix.rds"))
mat@threads <- threads

timing_cell <- system.time({
    cell_stats <- matrix_stats(mat, col_stats="variance")$col_stats
})

timing_gene <- system.time({
    gene_stats <- matrix_stats(mat, row_stats="variance")$row_stats
})

results <- tibble(
    time_cpu = c(sum(timing_cell[-3]), sum(timing_gene[-3])),
    time_elapsed = c(timing_cell[3], timing_gene[3]),
    axis = c("cell", "gene"),
)
write_tsv(results, output_timing)

write_tsv(
    tibble::as_tibble(t(cell_stats)),
    file.path(output_dir, "cell_stats.tsv")
)

write_tsv(
    tibble::as_tibble(t(gene_stats)),
    file.path(output_dir, "gene_stats.tsv")
)