suppressPackageStartupMessages({
    library(readr)
    library(tibble)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
reference_dir <- args[2]
staging_dir <- args[3]
output_timing <- args[4]

timing <- system.time({
    total_fragments <- system(sprintf("gunzip -c %s | wc -l", file.path(staging_dir, "sample.fragments.tsv.gz")), intern=TRUE)
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    total_fragments = total_fragments
)
write_tsv(results, output_timing)