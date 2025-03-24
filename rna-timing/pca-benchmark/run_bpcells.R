suppressPackageStartupMessages({
    library(BPCells)
    library(readr)
    library(tibble)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 6)

mode <- args[1]
sample_dir <- args[2]
data_dir <- args[3]
results_dir <- args[4]
tmp_dir <- args[5]
step <- args[6]

input_matrix <- file.path(sample_dir, "bpcells")
input_var_genes <- file.path(sample_dir, "variable_genes.txt")

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))

stopifnot(step %in% c("normalize", "pca"))
stopifnot(mode %in% c("stage_none", "stage_int", "stage_float", "randomized"))


normalized_mat <- file.path(tmp_dir, "normalize-mat")
normalized_rds <- file.path(tmp_dir, "normalize-mat.rds")
pca_rds <- file.path(data_dir, "pca.rds")

normalize_timing <- file.path(results_dir, "normalize-timing.tsv")
pca_timing <- file.path(results_dir, "pca-timing.tsv")

if (step == "normalize") {

    input_tmp <- file.path(tmp_dir, "input")
    system(sprintf("cp -r %s %s", input_matrix, input_tmp))
    
    mat <- open_matrix_dir(input_tmp)

    variable_genes <- readLines(input_var_genes)
    stopifnot(sum(rownames(mat) %in% variable_genes) == length(variable_genes))
    variable_genes <- rownames(mat)[rownames(mat) %in% variable_genes]

    
    
    if (mode == "stage_none" || mode == "randomized") {
        # Use the original data for everything
        timing <- system.time({
            stats <- matrix_stats(mat, row_stats="mean", col_stats="mean", threads=threads)
            row_sums <- stats$row_stats["mean",]*ncol(mat)
            col_sums <- stats$col_stats["mean",]*nrow(mat) 
            mat <- mat[row_sums > 0,]
            mat <- multiply_cols(mat, 10000/col_sums)
            
            # Just calculate mean and variance before taking the pre-calculated variable genes
            invisible(matrix_stats(mat, row_stats="variance", threads=threads))
            mat <- log1p(mat[variable_genes,])

            gene_stats <- matrix_stats(mat, row_stats="variance", threads=threads)$row_stats
            mat <- (mat - gene_stats["mean",]) / sqrt(gene_stats["variance",])
        })
        file_size <- tibble::tibble(
            file_bytes=sum(file.size(list.files(input_tmp, full.names=TRUE))),
            subset_nonzeros=sum(stats$row_stats["nonzero",variable_genes])
        )
    } else if(mode == "stage_float") {
        # Store the subset sparse normalized matrix to disk
        timing <- system.time({
            stats <- matrix_stats(mat, row_stats="mean", col_stats="mean", threads=threads)
            row_sums <- stats$row_stats["mean",]*ncol(mat)
            col_sums <- stats$col_stats["mean",]*nrow(mat) 
            mat <- mat[row_sums > 0,]
            mat <- multiply_cols(mat, 10000/col_sums)
            # Just calculate mean and variance before taking the pre-calculated variable genes
            invisible(matrix_stats(mat, row_stats="variance", threads=threads))
            mat <- log1p(mat[variable_genes,])

            gene_stats <- matrix_stats(mat, row_stats="variance", threads=threads)$row_stats
            mat <- mat / sqrt(gene_stats["variance",]) 

            # Save in floating-point format without dimnames to save some disk space and subsequent read bandwidth
            cellnames <- colnames(mat)
            colnames(mat) <- NULL
            mat <- mat |>
                convert_matrix_type("float") |>
                write_matrix_dir(normalized_mat)
            colnames(mat) <- cellnames
            mat <- mat - gene_stats["mean",] / sqrt(gene_stats["variance",])
        })
        file_size <- tibble::tibble(
            file_bytes=sum(file.size(list.files(normalized_mat, full.names=TRUE))),
            subset_nonzeros=sum(stats$row_stats["nonzero",variable_genes])
        )
    } else if(mode == "stage_int") {
        # Store the subset integer matrix to disk
        timing <- system.time({
            mat_raw <- mat
            stats <- matrix_stats(mat, row_stats="mean", col_stats="mean", threads=threads)
            row_sums <- stats$row_stats["mean",]*ncol(mat)
            col_sums <- stats$col_stats["mean",]*nrow(mat) 
            mat <- mat[row_sums > 0,]
            mat <- multiply_cols(mat, 10000/col_sums)
            # Just calculate mean and variance before taking the pre-calculated variable genes
            invisible(matrix_stats(mat, row_stats="variance", threads=threads))

            # Save without dimnames since this is just a temporary copy
            cellnames <- colnames(mat)
            colnames(mat_raw) <- NULL
            mat <- write_matrix_dir(mat_raw[variable_genes,], normalized_mat)
            colnames(mat) <- cellnames
            mat <- multiply_cols(mat, 10000/col_sums)
            mat <- log1p(mat)

            gene_stats <- matrix_stats(mat, row_stats="variance", threads=threads)$row_stats
            mat <- (mat - gene_stats["mean",]) / sqrt(gene_stats["variance",])
        })
        file_size <- tibble::tibble(
            file_bytes=sum(file.size(list.files(normalized_mat, full.names=TRUE))),
            subset_nonzeros=sum(stats$row_stats["nonzero",variable_genes])
        )
    } else {
        stop(paste0("Unrecognized mode: ", mode))
    }
    
    saveRDS(mat, normalized_rds, compress=FALSE)
    results <- tibble(
        time_cpu = sum(timing[-3]),
        time_elapsed = timing[3],
        step = step,
        n_ops = NA
    )
    write_tsv(results, normalize_timing)
    write_tsv(file_size, file.path(results_dir, "staged-file-size.tsv"))

} else if (step == "pca") {
    mat <- readRDS(normalized_rds)

    if (mode != "randomized") {
        timing <- system.time({
            svd <- svds(mat, k=50, threads=threads)
        })
    } else {
        source("../../utils/randomized_svd.R")
        timing <- system.time({
            mat <- BPCells:::parallel_split(mat, threads, threads*4)
            svd <- rsvd(mat, k=50, q=0)
        })
        svd$nops <- NA

    }

    saveRDS(svd, pca_rds, compress=FALSE)
    results <- tibble(
        time_cpu = sum(timing[-3]),
        time_elapsed = timing[3],
        step = step,
        n_ops = svd$nops
    )
    write_tsv(results, pca_timing)
}
