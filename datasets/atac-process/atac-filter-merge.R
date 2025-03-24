suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

project_dir <- args[1]
threads <- as.numeric(args[2])

dataset <- basename(project_dir)
ref_dir <- tempdir()

fragment_files <- Sys.glob(file.path(project_dir, "fragments", "*.fragments.tsv.gz"))
samples <- basename(str_remove(fragment_files, ".fragments.tsv.gz"))
fragment_dirs <- file.path(project_dir, "bpcells", samples)
prefix_barcodes <- "true" == read_lines(file.path(project_dir, "cell-barcodes-prefix.txt"))
genome <- read_lines(file.path(project_dir, "genome.txt"))
keeper_cells <- read_lines(file.path(project_dir, "cell-barcodes.txt"))

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))

open_fragments <- function(i, filter) {
    frags <- open_fragments_dir(fragment_dirs[i]) %>%
        select_chromosomes(standard_chrs)    
    if (prefix_barcodes) {
        frags <- prefix_cell_names(frags, paste0(samples[i], "."))
    }
    if (filter) {
        frags <- select_cells(frags, intersect(keeper_cells, cellNames(frags)))
    }
    return(frags)
}

write_filtered <- function(i) {
    open_fragments(i, filter=TRUE) %>%
        write_fragments_dir(file.path(project_dir, "bpcells_filtered", samples[i]))
}

cat(sprintf("Write filtered %s\n", Sys.time()))
fragment_list <- parallel::mclapply(
    mc.cores = threads,
    seq_along(samples),
    write_filtered
)

if (length(samples) > 1) {
    merge_input <- lapply(seq_along(samples), open_fragments, filter=FALSE) %>%
        do.call(c, .)
} else {
    merge_input <- open_fragments(1, filter=FALSE)
}

cat(sprintf("Write merge %s\n", Sys.time()))
merged_fragments <- write_fragments_dir(merge_input, file.path(project_dir, "bpcells", "merged"))

cat(sprintf("Write merge filtered %s\n", Sys.time()))
merged_fragments_filtered <- select_cells(merged_fragments, keeper_cells) %>%
    write_fragments_dir(file.path(project_dir, "bpcells_filtered", "merged"))

cat(sprintf("Success! %s\n", Sys.time()))