suppressPackageStartupMessages({
    library(BPCells)
    library(readr)
    library(tibble)
})


args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

tmp_dir <- args[1]
output_timing <- args[2]

mat <- open_matrix_dir(file.path(tmp_dir, "input"))

timing <- system.time({
    mat2 <- transpose_storage_order(mat, outdir=file.path(tmp_dir, "transpose"), tmpdir=tmp_dir)
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3]
)
write_tsv(results, output_timing)
