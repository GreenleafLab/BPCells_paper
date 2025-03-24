suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

input_type <- args[1]
census_root <- args[2]
output_dir <- args[3]
is_laptop <- as.logical(args[4])

stopifnot(input_type %in% c(
    "bpcells_cellmajor",
    "bpcells_genemajor",
    "bpcells_cellmajor_norm",
    "bpcells_genemajor_norm",
    "tiledb",
    "tiledb_norm"
))

if (str_detect(input_type, "bpcells_cellmajor")) {
    input_dir <- file.path(census_root, "bpcells", "subset_cellmajor")
} else if (str_detect(input_type, "bpcells_genemajor")) {
    input_dir <- file.path(census_root, "bpcells", "subset_genemajor")
} else if (input_type == "tiledb") {
    input_dir <- file.path(census_root, "tiledb", "raw_zstd_9")
} else if (input_type == "tiledb_norm") {
    input_dir <- file.path(census_root, "tiledb", "normalized_zstd_9")
} else {
    stop("Unhandled input_dir case")
}

if (is_laptop) {
    system(sprintf("ln -s %s %s", input_dir, file.path(output_dir, "matrix")))
} else {
    system(sprintf("cp -r %s %s", input_dir, file.path(output_dir, "matrix")))
}

if (str_detect(input_type, "bpcells")) {
    n <- length(list.files(file.path(output_dir, "matrix"), recursive=FALSE))
    paths <- file.path(file.path(output_dir, "matrix"), sprintf("%03d", seq_len(n)))
    mat_list <- lapply(paths, open_matrix_dir) |> lapply(convert_matrix_type, "double")

    mat <- do.call(cbind, mat_list)
    
    if (str_detect(input_type, "_norm")) {
        source("../normalize_gene_lengths.R")
        full_gene_mask <- readRDS(file.path(census_root, "obs_full_gene_mask.rds"))
        feature_length <- readRDS(file.path(census_root, "var.rds"))$feature_length
        cell_sums <- readRDS(file.path(census_root, "bpcells", "subset_norm", "raw_sum.rds"))
        mat <- normalize_gene_lengths(mat, full_gene_mask, feature_length) 
            
        # Do some awkwardness because BPCells doesn't yet automatically re-order the multiply_cols operation with earlier merges.
        # Therefore, we need to apply multiply_cols prior to merging in order to get parallel reads when the inputs are manually concatenated
        cell_counts <- vapply(mat@matrix_list, ncol, integer(1))
        start_idx <- cumsum(c(0,cell_counts[-length(cell_counts)])) + 1
        mat@matrix_list <- lapply(seq_along(mat@matrix_list), function(i) {
            multiply_cols(mat@matrix_list[[i]], 1/cell_sums[start_idx[i]:(start_idx[i]+cell_counts[i]-1)])
        })
    } 
    saveRDS(mat, file.path(output_dir, "matrix.rds"), compress=FALSE)
}