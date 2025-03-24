suppressPackageStartupMessages({
    library(tidyverse)
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

dir <- args[1]
replicate <- args[2]
output_path <- args[3]

m <- open_matrix_dir(dir) |>
    write_matrix_memory(compress=FALSE)
mc <- write_matrix_memory(m)

n <- length(m@val)

val_data <- rep_len(12321L, n)
val_idx <- rep_len(12321L, n/64)
val_idx_offsets <- rep_len(1.123, n/64)
val_result <- rep_len(12321L, n + 128)

index_data <- rep_len(12321L, n)
index_idx <- rep_len(12321L, n/64)
index_idx_offsets <- rep_len(1.123, n/64)
index_starts <- rep_len(12321L, n/64)
index_result <- rep_len(12321L, n + 128)

# Test for round-trip correctness, then iterate + time
# Note: idx_offsets is actually stored slightly different because we don't cast back to a double
equal_prefix <- function(prefix, data) {
    all.equal(prefix, data[seq_along(prefix)])
}

BPCells:::write_bp128_for(m@val, val_data, val_idx, val_idx_offsets)
BPCells:::read_bp128_for(val_data, val_idx, val_idx_offsets, val_result, n)
stopifnot(equal_prefix(m@val, val_result))
stopifnot(equal_prefix(mc@val_data, val_data))
stopifnot(equal_prefix(mc@val_idx, val_idx))

BPCells:::write_bp128_d1z(m@index, index_data, index_idx, index_idx_offsets, index_starts)
BPCells:::read_bp128_d1z(index_data, index_idx, index_idx_offsets, index_starts, index_result, n)
stopifnot(equal_prefix(m@index, index_result))
stopifnot(equal_prefix(mc@index_data, index_data))
stopifnot(equal_prefix(mc@index_idx, index_idx))
stopifnot(equal_prefix(mc@index_starts, index_starts))


bench::mark(
    write_val = BPCells:::write_bp128_for(m@val, val_data, val_idx, val_idx_offsets),
    read_val = BPCells:::read_bp128_for(val_data, val_idx, val_idx_offsets, val_result, n),
    write_index = BPCells:::write_bp128_d1z(m@index, index_data, index_idx, index_idx_offsets, index_starts),
    read_index = BPCells:::read_bp128_d1z(index_data, index_idx, index_idx_offsets, index_starts, index_result, n),
    min_time = 2
) |>
    select(!c(result, memory, time, gc)) |>
    mutate(
        min=as.numeric(min),
        median=as.numeric(median),
        input = dir, 
        replicate=replicate
    ) |>
    write_tsv(output_path)

