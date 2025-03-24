suppressPackageStartupMessages({
    library(tidyverse)
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

dir <- args[1]
replicate <- args[2]
output_path <- args[3]

f <- open_fragments_dir(dir) |>
    write_fragments_memory(compress=FALSE)
fc <- write_fragments_memory(f)

n <- length(f@cell)

cell_data <- rep_len(12321L, n)
cell_idx <- rep_len(12321L, n/64)
cell_idx_offsets <- rep_len(1.123, n/64)
cell_result <- rep_len(12321L, n + 128)

start_data <- rep_len(12321L, n)
start_idx <- rep_len(12321L, n/64)
start_idx_offsets <- rep_len(1.123, n/64)
start_starts <- rep_len(12321L, n/64)
start_result <- rep_len(12321L, n + 128)

end_data <- rep_len(12321L, n)
end_idx <- rep_len(12321L, n/64)
end_idx_offsets <- rep_len(1.123, n/64)
end_result <- rep_len(12321L, n + 128)

# Test for round-trip correctness, then iterate + time
# Note: idx_offsets is actually stored slightly different because we don't cast back to a double
equal_prefix <- function(prefix, data) {
    all.equal(prefix, data[seq_along(prefix)])
}

BPCells:::write_bp128(f@cell, cell_data, cell_idx, cell_idx_offsets)
BPCells:::read_bp128(cell_data, cell_idx, cell_idx_offsets, cell_result, n)
stopifnot(equal_prefix(f@cell, cell_result))
stopifnot(equal_prefix(fc@cell_data, cell_data))
stopifnot(equal_prefix(fc@cell_idx, cell_idx))

BPCells:::write_bp128_d1(f@start, start_data, start_idx, start_idx_offsets, start_starts)
BPCells:::read_bp128_d1(start_data, start_idx, start_idx_offsets, start_starts, start_result, n)
stopifnot(equal_prefix(f@start, start_result))
stopifnot(equal_prefix(fc@start_data, start_data))
stopifnot(equal_prefix(fc@start_idx, start_idx))
stopifnot(equal_prefix(fc@start_starts, start_starts))

BPCells:::write_bp128_end(f@end, f@start, end_data, end_idx, end_idx_offsets)
BPCells:::read_bp128_end(end_data, end_idx, end_idx_offsets, f@start, end_result, n)
stopifnot(equal_prefix(f@end, end_result))
stopifnot(equal_prefix(fc@end_data, end_data))
stopifnot(equal_prefix(fc@end_idx, end_idx))

bench::mark(
    write_cell = BPCells:::write_bp128(f@cell, cell_data, cell_idx, cell_idx_offsets),
    read_cell = BPCells:::read_bp128(cell_data, cell_idx, cell_idx_offsets, cell_result, n),
    write_start = BPCells:::write_bp128_d1(f@start, start_data, start_idx, start_idx_offsets, start_starts),
    read_start = BPCells:::read_bp128_d1(start_data, start_idx, start_idx_offsets, start_starts, start_result, n),
    write_end = BPCells:::write_bp128_end(f@end, f@start, end_data, end_idx, end_idx_offsets),
    read_end = BPCells:::read_bp128_end(end_data, end_idx, end_idx_offsets, f@start, end_result, n),
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

