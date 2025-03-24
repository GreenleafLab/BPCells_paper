suppressPackageStartupMessages({
    library(tiledbsoma)
    library(tidyverse)
    library(BPCells)
    library(Matrix)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

tiledb_path <- args[1]
output_base <- args[2]
chunk_num <- as.integer(args[3])

x <- SOMASparseNDArrayOpen(file.path(tiledb_path, "ms/RNA/X/raw"))

# Rows are cells
cells <- x$shape()[1]
genes <- x$shape()[2]

r <- x$tiledb_array()

read_matrix_chunk <- function(mat, row_start, row_count) {
    rows <- as.integer(row_start + 0:(row_count-1))
    rows <- rows[rows < cells]

    data <- mat[rows,]
    if (length(data$soma_dim_1) > 1.5e9) {
        msg <- sprintf("Read chunk size too large (%d) from row (%d) chunk size (%d)", length(data$soma_dim_1), row_start, row_count)
        cat(msg)
        stop(msg)
    }
    if (any(floor(data$soma_data) != data$soma_data)) {
        msg <- sprintf("Found non-integer data from row (%d) chunk size (%d)", row_start, row_count)
        cat(msg)
        stop(msg)
    }
    sparseMatrix(
        i = data$soma_dim_1 + 1L,
        j = data$soma_dim_0 - min(rows) + 1L,
        x = data$soma_data,
        dims = c(genes, length(rows))
    ) %>%
        as("IterableMatrix") %>%
        convert_matrix_type("uint32_t") %>%
        write_matrix_memory()
}

# Read 10k cells at a time into a dgCMatrix, store in-memory compressed, then concatenate into 100K cell chunks for saving
export_matrix_chunk <- function(chunk_num, mat, out_base_dir) {
    chunk_size <- 100000
    batches <- 10
    batch_size <- chunk_size/batches
    stopifnot(floor(batch_size) == batch_size)
    
    mat_batches <- list()
    for (i in seq_len(batches)) {
        row_start <- chunk_size * (chunk_num - 1) + batch_size * (i-1)
        if (row_start >= cells) {
            break
        }
        batch_mat <- read_matrix_chunk(mat, row_start=row_start, row_count = batch_size)
        mat_batches <- c(mat_batches, batch_mat)
        gc()
    }
    
    out_path_cell <- file.path(out_base_dir, sprintf("cellmajor_chunks/%03d", chunk_num))
    do.call(cbind, mat_batches) %>%
        write_matrix_dir(out_path_cell)

    return(NULL)
}

system.time({
    export_matrix_chunk(chunk_num, r, output_base)
})
