suppressPackageStartupMessages({
    library(BPCells)
    library(Matrix)
    library(tidyverse)
})

script_dir <- commandArgs(trailingOnly=FALSE)
script_dir <- script_dir[str_detect(script_dir, "--file=")] |> str_replace("--file=", "") |> dirname() |> normalizePath()
source(file.path(script_dir, "../bitwidth-helpers.R"))

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)
input_dir <- args[1]
tmp_dir <- args[2]
output_dir <- args[3]

mat <- open_matrix_dir(input_dir)
mat_orderings <- list(
    "original" = input_dir,
    "mean" = file.path(tmp_dir, "matrix_sortmean"),
    "nonzero" = file.path(tmp_dir, "matrix_sortnnz"),
    "shuffle" = file.path(tmp_dir, "matrix_shuffle")
)



is_colmajor <- storage_order(mat) == "col"
if (is_colmajor) {
    stats <- matrix_stats(mat, row_stats="mean")$row_stats
} else {
    stats <- matrix_stats(mat, col_stats="mean")$col_stats
    mat <- t(mat)
}
matrix_sortmean <- mat[order(stats["mean",]), ]
matrix_sortnnz <- mat[order(stats["nonzero",]), ] 
set.seed(125124)
matrix_shuffle <- mat[sample.int(nrow(mat)), ]
rm(mat)
gc()
if (!is_colmajor) {
    matrix_sortmean <- t(matrix_sortmean)
    matrix_sortnnz <- t(matrix_sortnnz)
    matrix_shuffle <- t(matrix_shuffle)
}
write_matrix_dir(matrix_sortmean, mat_orderings[["mean"]])
rm(matrix_sortmean)
gc()
write_matrix_dir(matrix_sortnnz, mat_orderings[["nonzero"]])
rm(matrix_sortnnz)
gc()
write_matrix_dir(matrix_shuffle, mat_orderings[["shuffle"]])
rm(matrix_shuffle)
gc()


# Calculate matrix bitwidths
bitwidths <- lapply(mat_orderings, calculate_bitwidths_matrix) |>
    bind_rows(.id = "ordering") |>
    select(!folder)


write_tsv(bitwidths, file.path(output_dir, "bitwidths.tsv"))

# Calculate storage usage
storage <- lapply(mat_orderings, matrix_storage) |>
    bind_rows(.id = "ordering") |>
    select(!folder)


write_tsv(storage, file.path(output_dir, "storage.tsv"))

# Calculate value percentiles
value_histogram <- matrix_value_histogram(input_dir)

write_tsv(value_histogram, file.path(output_dir, "value_histogram.tsv"))

