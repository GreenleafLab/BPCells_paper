suppressPackageStartupMessages({
    library(tidyverse)
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_path <- args[1]
output_path <- args[2]

open_fragments_10x(input_path) %>%
    write_fragments_dir(output_path)