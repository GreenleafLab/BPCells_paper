# After convert_bpcells_chunks.R, make subsets
suppressPackageStartupMessages({
    library(BPCells)
    library(readr)
    library(tibble)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 5)

input_dir <- args[1]
census_root <- args[2]
total_chunks <- as.integer(args[3])
output_dir <- args[4]
output_timing <- args[5]

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
obs_idx <- readRDS(file.path(census_root, "obs_idx.rds"))

system.time({
    # Check we aren't missing any input chunks
    input_chunks <- as.integer(list.dirs(input_dir, recursive=FALSE, full.names=FALSE))
    stopifnot(all(input_chunks == seq_along(input_chunks)))

    # Load the matrix
    mat_list <- list.files(input_dir, include.dirs=TRUE, recursive=FALSE, full.names=TRUE) |>
        lapply(open_matrix_dir)
    mat <- do.call(cbind, mat_list)
})

write_chunk <- function(chunk_num) {
    # Subset in chunks of 100k
    chunk_size <- 100000
    rows <- 0:(chunk_size-1)
    rows <- rows + chunk_size * (chunk_num - 1)
    rows <- rows[rows < length(obs_idx)]
    rows <- as.integer(rows + 1)

    # Subset to just the cells we want in this chunk
    mat[, obs_idx[rows]] |>
        write_matrix_dir(sprintf("%s/%03d", output_dir, chunk_num))
}

timing <- system.time({
    parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE,
        seq_len(total_chunks),
        write_chunk
    )
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    bytes = sum(file.size(list.files(output_dir, recursive=TRUE, full.names=TRUE)))
)

write_tsv(results, output_timing)
