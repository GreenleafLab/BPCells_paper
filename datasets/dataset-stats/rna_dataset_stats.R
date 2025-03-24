suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_matrix <- args[1]
output_stats <- args[2]

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))

mat <- open_matrix_dir(input_matrix) 
stats <-  matrix_stats(mat, col_stats="mean", threads=threads) 

reads <- stats$col_stats["mean",] * nrow(mat)

# dataset	cells	genes	median_reads	median_genes	nonzero_entries	fraction_nonzero	mean_nonzero_val
total_nonzero <- sum(stats$col_stats["nonzero",])

tibble::tibble(
    cells = ncol(mat),
    genes = nrow(mat),
    median_reads = median(reads),
    median_genes = median(stats$col_stats["nonzero",]),
    nonzero_entries = total_nonzero,
    fraction_nonzero = as.numeric(total_nonzero) / as.numeric(cells) / as.numeric(genes),
    mean_nonzero_val = sum(reads)/total_nonzero
) |> write_tsv(output_stats)
