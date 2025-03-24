suppressPackageStartupMessages({
    library(BPCells)
    library(tibble)
    library(readr)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 5)

input_dir <- args[1]
coords <- args[2]
axis <- args[3]
output_timing <- args[4]
threads <- as.integer(args[5])

stopifnot(axis %in% c("cell", "gene"))

mat <- readRDS(file.path(input_dir, "matrix.rds"))
mat@threads <- threads
coords <- readLines(coords) |> as.integer()

# Count nonzeros along the most major storage axis
if (mat@transpose) {
    count_nonzeros <- function(m) {
        sum(matrix_stats(m, row_stats="nonzero")$row_stats)
    }
} else {
    count_nonzeros <- function(m) {
        sum(matrix_stats(m, col_stats="nonzero")$col_stats)
    }
}

if (axis == "gene" && length(coords) < nrow(mat)) {
    cat("Subsetting on gene axis\n")
    timing <- system.time({
        total_nonzero <- count_nonzeros(mat[coords,])
    })
} else if (axis == "cell" && length(coords) < ncol(mat)) {
    cat("Subsetting on cell axis\n")
    timing <- system.time({
        total_nonzero <- count_nonzeros(mat[,coords])
    })
} else {
    cat("No subsetting\n")
    timing <- system.time({
        total_nonzero <- count_nonzeros(mat)
    })
}

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    entries_loaded = total_nonzero
)
write_tsv(results, output_timing)