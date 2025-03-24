suppressPackageStartupMessages({
    library(readr)
    library(tibble)
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
reference_dir <- args[2]
staging_dir <- args[3]
output_timing <- args[4]

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
timing <- system.time({
    frags <- open_fragments_10x(file.path(staging_dir, "sample.fragments.tsv.gz")) |>
        select_chromosomes(standard_chrs) |>
        write_fragments_dir(file.path(staging_dir, "sample"))
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    bytes = sum(file.size(list.files(file.path(staging_dir, "sample"), full.names=TRUE)))
)

write_tsv(results, output_timing)