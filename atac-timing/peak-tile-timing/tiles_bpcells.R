suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
output_timing <- args[2]
staging_dir <- args[3]
ref_dir <- args[4]

genome <- readLines(file.path(dataset_dir, "genome.txt"))

if (genome == "hg19") {
    chr_sizes <- read_ucsc_chrom_sizes(ref_dir, "hg19")
} else if (genome == "hg38") {
    chr_sizes <- read_ucsc_chrom_sizes(ref_dir, "hg38")
}

chrs <- open_fragments_dir(file.path(staging_dir, "sample")) |> chrNames()
tiles <- chr_sizes %>%
    mutate(tile_width=500L) %>%
    filter(chr %in% chrs)
tiles <- tiles[order_ranges(tiles, chrs),]

cat(sprintf("Calculating matrices %s\n", Sys.time()))

timing <- system.time({
    open_fragments_dir(file.path(staging_dir, "sample")) |>
        tile_matrix(tiles) |>
        write_matrix_dir(file.path(staging_dir, "matrix"))
})

cat(sprintf("Outputting timing results %s\n", Sys.time()))

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    tool = "bpcells",
)

write_tsv(results, output_timing)
cat(sprintf("Done %s\n", Sys.time()))
