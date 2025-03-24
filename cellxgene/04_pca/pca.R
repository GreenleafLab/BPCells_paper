suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 7)

input_mat <- args[1]
census_root <- args[2]
work_dir <- args[3]
compress <- args[4]
precision <- args[5]
time_output <- args[6]
pca_output <- args[7]

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))

stopifnot(compress %in% c("compress", "uncompress"))
stopifnot(precision %in% c("exact", "randomized"))

n <- length(list.files(input_mat, recursive=FALSE))
mat_list <- parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, seq_len(n), function(i) {
    return(open_matrix_dir(file.path(input_mat, sprintf("%03d", i))) |> convert_matrix_type("double"))
})

mat <- do.call(cbind, mat_list)
mat@threads <- threads

chunk_size <- 100000
total_cells <- ncol(mat)
total_chunks <- ceiling(total_cells / chunk_size)

chunk_indices <- function(chunk) {
    stopifnot(ncol(mat_list[[chunk]]) == chunk_size || chunk == total_chunks)
    rows <- 0:(chunk_size-1)
    rows <- rows + chunk_size * (chunk - 1)
    rows <- rows[rows < ncol(mat)]
    as.integer(1 + rows)
}

# Perform normalization
# 
# The unusual structure of this code is because BPCells does not currently
# implement operations to postpone merge operations until after transformation operations.
#
# In order to achieve parallelism in cases like this where we are trying to explicitly
# parallelize by pre-determined file boundaries, then this approach is the current best option,
# though future improvements will hopefully enable a more natural coding style
cat(sprintf("Normalization part 1 %s\n", Sys.time()))
source("../normalize_gene_lengths.R")
full_gene_mask <- readRDS(file.path(census_root, "obs_full_gene_mask.rds"))
feature_length <- readRDS(file.path(census_root, "var.rds"))$feature_length
system.time({
    mat <- normalize_gene_lengths(mat, full_gene_mask, feature_length)
    cell_counts <- vapply(mat@matrix_list, ncol, integer(1))
    start_idx <- cumsum(c(0, cell_counts[-length(cell_counts)])) + 1
    end_idx <- cumsum(cell_counts)

    stats <- matrix_stats(mat, col_stats="mean")

    # mat <- multiply_cols(mat, 10000/(nrow(mat)*stats$col_stats["mean",])) |>
    #     log1p()
    mat@matrix_list <- parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, seq_along(mat@matrix_list), function(i) {
        col_sums <- stats$col_stats["mean", start_idx[i]:end_idx[i]] * nrow(mat)
        m <- mat@matrix_list[[i]] |>
            multiply_cols(10000/col_sums) |>
            log1p()
        m
    })

    gene_stats <- matrix_stats(mat, row_stats="variance")$row_stats
    
    gene_sd <- sqrt(gene_stats["variance",])

    # Put in a dummy value so we can handle all-0 genes smoothly (seem to cover 232 genes in our matrix)
    gene_sd[gene_sd == 0] <- 0.01
})

if (compress == "uncompress") {
    cat(sprintf("Saving normalized matrix to disk %s\n", Sys.time()))
    mat@matrix_list <- parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, seq_along(mat@matrix_list), function(i) {
        m <- convert_matrix_type(mat@matrix_list[[i]], "float")
        write_matrix_dir(m, file.path(work_dir, "uncompress_mat", sprintf("%03d", i)), compress=FALSE)
    })
}

cat(sprintf("Normalization part 2 %s\n", Sys.time()))
system.time({
    mat@matrix_list <- parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, seq_along(mat@matrix_list), function(i) {
        m <- mat@matrix_list[[i]]
        m <- (m - gene_stats["mean",]) / gene_sd  
        m  
    })
})


# Perform PCA
cat(sprintf("PCA %s\n", Sys.time()))

if (precision == "exact") {
    timing <- system.time({
        svd_res <- svds(mat, k=32)
    })
} else {
    source("../../utils/randomized_svd.R")
    timing <- system.time({
        svd_res <- rsvd(mat, k=32)
    })
}

cat(sprintf("Saving results %s\n", Sys.time()))

saveRDS(svd_res, pca_output, compress=FALSE)

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    n_ops=svd_res$nops
)
write_tsv(results, time_output)

cat(sprintf("Done %s\n", Sys.time()))
