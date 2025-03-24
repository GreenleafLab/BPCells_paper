suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_dir <- args[1]
output_timing <- args[2]

sample_paths <- list.files(input_dir, recursive=FALSE, full.names=TRUE)

timing <- system.time({
    merge_input <- lapply(sample_paths, open_fragments_dir) %>% do.call(c, .)
    merged_fragments <- write_fragments_dir(merge_input, file.path(input_dir, "merged"))
})


results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3]
)
write_tsv(results, output_timing)