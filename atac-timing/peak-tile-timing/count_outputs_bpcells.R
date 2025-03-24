suppressPackageStartupMessages({
    library(readr)
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

staging_dir <- args[1]
matrix_type <- args[2]
output_path <- args[3]

stopifnot(matrix_type %in% c("peak", "tile"))

cat(sprintf("Counting overlaps %s\n", Sys.time()))

mat <- open_matrix_dir(file.path(staging_dir, "matrix"))
total_overlaps <- sum(rowSums(mat))

write_lines(total_overlaps, output_path)

cat(sprintf("Done %s\n", Sys.time()))