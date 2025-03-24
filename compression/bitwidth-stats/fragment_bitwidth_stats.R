suppressPackageStartupMessages({
    library(BPCells)
    library(Matrix)
    library(tidyverse)
})

script_dir <- commandArgs(trailingOnly=FALSE)
script_dir <- script_dir[str_detect(script_dir, "--file=")] |> str_replace("--file=", "") |> dirname() |> normalizePath()
source(file.path(script_dir, "../bitwidth-helpers.R"))

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)
input_dir <- args[1]
output_dir <- args[2]

samples <- list.dirs(input_dir, recursive=FALSE, full.names=TRUE)
names(samples) <- basename(samples)

# Calculate fragment bitwidths
bitwidths <- lapply(samples, calculate_bitwidths_fragments) |>
    bind_rows(.id = "sample") |>
    select(!folder)

write_tsv(bitwidths, file.path(output_dir, "bitwidths.tsv"))

# Calculate storage usage
storage <- lapply(samples, fragments_storage) |>
    bind_rows(.id = "sample") |>
    select(!folder)

write_tsv(storage, file.path(output_dir, "storage.tsv"))