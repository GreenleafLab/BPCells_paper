suppressPackageStartupMessages({
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

input_path <- args[1]
file_type <- args[2]
output_path <- args[3]

open_dir <- if(file_type == "fragments") open_fragments_dir else open_matrix_dir
write_dir <- if(file_type == "fragments") write_fragments_dir else write_matrix_dir

# Prep an uncompressed version of a dataset
open_dir(input_path) |> write_dir(output_path, compress=FALSE)