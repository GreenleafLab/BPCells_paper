suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

project_dir <- args[1]
staging_dir <- args[2]

dataset <- basename(project_dir)
fragment_files <- Sys.glob(file.path(project_dir, "fragments", "*.fragments.tsv.gz"))
samples <- str_remove(basename(fragment_files), ".fragments.tsv.gz")
fragment_dirs <- file.path(project_dir, "bpcells_filtered", samples)

file.copy(fragment_dirs, staging_dir, recursive=TRUE)
