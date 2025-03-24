suppressPackageStartupMessages({
    library(Matrix)
    library(readr)
    library(tibble)
})


args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

tmp_dir <- args[1]
output_timing <- args[2]

h5 <- hdf5r::H5File$new(file.path(tmp_dir, "input.h5"))
mat <- sparseMatrix(
    i =  h5[["matrix/indices"]][] + 1,
    p = h5[["matrix/indptr"]][],
    x = h5[["matrix/data"]][]
)

gc()

timing <- system.time({
    mat2 <- t(mat)
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3]
)
write_tsv(results, output_timing)
