suppressPackageStartupMessages({
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

input_peaks <- args[1]
output_peaks <- args[2]
num_peaks <- as.integer(args[3])
chr_names_file <- args[4]


peaks <- read_tsv(input_peaks, col_names=c("chr", "start", "end"), col_types="cii") 

valid_chromosomes <- read_lines(chr_names_file)
# Note that 3k pbmc only has like 98k peaks so we'll sometimes not be able to get the full 100k
set.seed(1251234)
peaks |>
    distinct() |>
    filter(chr %in% valid_chromosomes) |>
    slice_sample(n=min(num_peaks, nrow(peaks))) |>
    arrange(chr, start) |>
    write_tsv(output_peaks, col_names=FALSE)