suppressPackageStartupMessages({
    library(BPCells)
    library(readr)
    library(tibble)
})

source("../normalize_gene_lengths.R")

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 5)

input_dir <- args[1]
census_root <- args[2]
total_chunks <- as.integer(args[3])
output_dir <- args[4]
output_timing <- args[5]

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
obs_idx <- readRDS(file.path(census_root, "obs_idx.rds"))
full_gene_mask <- readRDS(file.path(census_root, "obs_full_gene_mask.rds"))
feature_length <- readRDS(file.path(census_root, "var.rds"))$feature_length

system.time({
    # Check we aren't missing any input chunks
    input_chunks <- as.integer(list.dirs(input_dir, recursive=FALSE, full.names=FALSE))
    stopifnot(all(input_chunks == seq_along(input_chunks)))

    # Load the matrix
    mat_list <- list.files(input_dir, include.dirs=TRUE, recursive=FALSE, full.names=TRUE) |>
        lapply(open_matrix_dir) |>
        lapply(convert_matrix_type, "double")
    mat <- do.call(cbind, mat_list)

    stopifnot(ncol(mat) > length(full_gene_mask))
    # Make an all-cell version of the full_gene_mask
    extended_mask <- rep.int(FALSE, ncol(mat))
    extended_mask[obs_idx] <- full_gene_mask
    
    mat@threads <- threads
    mat <- normalize_gene_lengths(mat, extended_mask, feature_length)
})

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
timing <- system.time({
    saveRDS(colSums(mat)[obs_idx], file.path(output_dir, "raw_sum.rds"), compress=FALSE)
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    bytes = sum(file.size(list.files(output_dir, recursive=TRUE, full.names=TRUE)))
)

write_tsv(results, output_timing)
