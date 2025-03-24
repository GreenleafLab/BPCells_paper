suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

project_dir <- args[1]
threads <- as.numeric(args[2])
ref_dir <- args[3]


fragment_files <- Sys.glob(file.path(project_dir, "fragments", "*.fragments.tsv.gz"))
samples <- basename(str_remove(fragment_files, ".fragments.tsv.gz"))
fragment_dirs <- file.path(project_dir, "bpcells", samples)
prefix_barcodes <- "true" == read_lines(file.path(project_dir, "cell-barcodes-prefix.txt"))
genome <- read_lines(file.path(project_dir, "genome.txt"))
keeper_cells <- read_lines(file.path(project_dir, "cell-barcodes.txt"))

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
peaks <- read_tsv(file.path(project_dir, "peaks.bed"), col_names=c("chr", "start", "end"), col_types="cii") %>%
    filter(chr %in% standard_chrs)

peaks <- peaks[order_ranges(peaks, standard_chrs),]

stopifnot(nrow(distinct(peaks)) == nrow(peaks))

if (genome == "hg19") {
    chr_sizes <- read_ucsc_chrom_sizes(ref_dir, "hg19")
} else if (genome == "hg38") {
    chr_sizes <- read_ucsc_chrom_sizes(ref_dir, "hg38")
}

tiles <- chr_sizes %>%
    mutate(tile_width=500L)
tiles <- tiles[order_ranges(tiles, standard_chrs),]

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

write_peaks <- function(i) {
    open_fragments(i, TRUE) %>%
        peak_matrix(peaks) %>%
        write_matrix_dir(file.path(project_dir, "matrices", "peaks", samples[i]))
}

write_tiles <- function(i) {
    res <- open_fragments(i, TRUE) %>%
        tile_matrix(tiles) %>%
        write_matrix_dir(file.path(project_dir, "matrices", "tiles", samples[i]))
    res@dimnames <- list(NULL, res@dimnames[[2]])
    res
}

cat(sprintf("Write peak mats %s\n", Sys.time()))
peak_mats <- parallel::mclapply(mc.cores=threads, seq_along(samples), write_peaks)
if (length(peak_mats) > 1) {
    peaks_merged <- do.call(cbind, peak_mats)
} else {
    peaks_merged <- peak_mats[[1]]
}
cat(sprintf("Write peaks merged %s\n", Sys.time()))
merge_peak_mat <- write_matrix_dir(peaks_merged, file.path(project_dir, "matrices", "peaks", "merged"))

# Avoid memory usage issues on many-sample datasets
rm(merge_peak_mat, peaks_merged, peak_mats)
gc()

cat(sprintf("Write tile mats %s\n", Sys.time()))
tile_mats <- parallel::mclapply(mc.cores=threads, seq_along(samples), write_tiles)
if (length(tile_mats) > 1) {
    tiles_merged <- do.call(cbind, tile_mats)
} else {
    tiles_merged <- tile_mats[[1]]
}
rownames(tiles_merged) <- rownames(open_matrix_dir(file.path(project_dir, "matrices", "tiles", samples[1])))
cat(sprintf("Write tiles merged %s\n", Sys.time()))
merge_tile_mat <- write_matrix_dir(tiles_merged, file.path(project_dir, "matrices", "tiles", "merged"))

# Avoid memory usage issues on many-sample datasets
rm(merge_tile_mat, tiles_merged, tile_mats)
gc()


cat(sprintf("Re-ordering peak matrices %s\n", Sys.time()))
parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, c(samples, "merged"), function (s) {
    open_matrix_dir(file.path(project_dir, "matrices", "peaks", s)) |>
        transpose_storage_order(file.path(project_dir, "matrices", "peaks_transpose", s))
});
gc()
cat(sprintf("Re-ordering tile matrices %s\n", Sys.time()))
parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, c(samples, "merged"), function (s) {
    open_matrix_dir(file.path(project_dir, "matrices", "tiles", s)) |>
        transpose_storage_order(file.path(project_dir, "matrices", "tiles_transpose", s))
});
gc()

cat(sprintf("Success! %s\n", Sys.time()))
