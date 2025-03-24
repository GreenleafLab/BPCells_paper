suppressPackageStartupMessages({
    library(readr)
    library(ArchR)
})


args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

staging_dir <- args[1]
matrix_type <- args[2]
output_path <- args[3]

stopifnot(matrix_type %in% c("peak", "tile"))

if (matrix_type == "peak") {
    matrix_name <- "PeakMatrix"
} else {
    matrix_name <- "TileMatrix"
}

proj <- readRDS(file.path(staging_dir, "archr_proj.rds"))
matrix_sum <- sum(ArchR:::.getRowSums(proj@sampleColData$ArrowFiles, matrix_name, threads=1L)$rowSums)

write_lines(matrix_sum, output_path)