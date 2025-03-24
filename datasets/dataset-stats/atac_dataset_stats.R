suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

project_dir <- args[1]
output_stats <- args[2]

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))

fragment_files <- Sys.glob(file.path(project_dir, "fragments", "*.fragments.tsv.gz"))
samples <- basename(str_remove(fragment_files, ".fragments.tsv.gz"))
fragment_dirs <- file.path(project_dir, "bpcells", samples)
prefix_barcodes <- "true" == read_lines(file.path(project_dir, "cell-barcodes-prefix.txt"))
genome <- read_lines(file.path(project_dir, "genome.txt"))
keeper_cells <- read_lines(file.path(project_dir, "cell-barcodes.txt"))

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))


# total_barcodes - pre-filtering
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

dataset_stats <- function(i, filter) {
    frags <- open_fragments(i, filter)
    frags_per_cell <- nucleosome_counts(frags, nucleosome_width=1e6)$nFrags 

    tibble::tibble(
        sample = samples[i],
        filtered = filter,
        barcodes=length(cellNames(frags)),
        total_fragments=sum(frags_per_cell),
        median_fragments=median(frags_per_cell),
        frags_per_cell=list(frags_per_cell)
    )
}

args <- tidyr::expand_grid(i=seq_along(samples), filtered=c(TRUE, FALSE))
results <- parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, seq_len(nrow(args)), function(i) dataset_stats(args$i[i], args$filtered[i])) |>
    dplyr::bind_rows()

merged_stats <- results |>
    group_by(filtered) |>
    summarize(
        sample="merged",
        barcodes=sum(barcodes),
        total_fragments=sum(total_fragments),
        median_fragments=median(unlist(frags_per_cell))
    ) |>
    select(sample, filtered, barcodes, total_fragments, median_fragments)

dplyr::bind_rows(select(results, !frags_per_cell), merged_stats) |>
    write_tsv(output_stats)