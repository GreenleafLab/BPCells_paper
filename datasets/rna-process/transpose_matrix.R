suppressPackageStartupMessages({
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_dir <- args[1]
output_dir <- args[2]

cat(sprintf("Transposing matrix start: %s\n", Sys.time()))
open_matrix_dir(input_dir) |>
    transpose_storage_order(outdir=output_dir)
cat(sprintf("Transposing matrix end: %s\n", Sys.time()))
