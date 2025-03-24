suppressPackageStartupMessages({
    library(BPCells)
    library(magrittr)
    library(Matrix)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_path <- args[1]
output_path <- args[2]

stopifnot(grepl(".h5$", input_path))

open_matrix_10x_hdf5(input_path) %>%
    convert_matrix_type("uint32_t") %>%
    write_matrix_dir(output_path);

