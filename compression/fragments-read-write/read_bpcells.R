# Import *only* fragments from ArchR, without running QC
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

timing <- system.time({
    # Note: a faster way to simply scan over fragments on disk is with nucleosome_counts(),
    # but here we load all data into a copy in memory since other tool benchmarks
    # have to fall back on that option.
    in_memory_fragments <- open_fragments_dir(file.path(staging_dir, "sample")) |>
        write_fragments_memory()
})

total_fragments <- nucleosome_counts(in_memory_fragments)$nFrags |> sum()

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    total_fragments = total_fragments
)
write_tsv(results, output_timing)