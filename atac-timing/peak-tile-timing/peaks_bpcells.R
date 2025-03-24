suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
output_timing <- args[2]
peaks_file <- args[3]
staging_dir <- args[4]

cat(sprintf("Loading peaks %s\n", Sys.time()))

chrs <- open_fragments_dir(file.path(staging_dir, "sample")) |> chrNames()
peaks <- read_tsv(peaks_file, col_names=c("chr", "start", "end"), col_types="cii") |>
    filter(chr %in% chrs)
peaks <- peaks[order_ranges(peaks, chrs),]
stopifnot(nrow(distinct(peaks)) == nrow(peaks))

cat(sprintf("Calculating matrices %s\n", Sys.time()))


timing <- system.time({
    open_fragments_dir(file.path(staging_dir, "sample")) |>
        peak_matrix(peaks) |>
        write_matrix_dir(file.path(staging_dir, "matrix"))
})

cat(sprintf("Outputting timing results %s\n", Sys.time()))
results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    tool = "bpcells",
    peaks = nrow(peaks)
)
write_tsv(results, output_timing)

cat(sprintf("Done %s\n", Sys.time()))
