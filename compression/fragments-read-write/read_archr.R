# Import *only* fragments from ArchR, without running QC
suppressPackageStartupMessages({
    library(readr)
    library(tibble)
    library(ArchR)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
reference_dir <- args[2]
staging_dir <- args[3]
output_timing <- args[4]

staging_input <- file.path(staging_dir, "sample.arrow")
total_fragments <- 0L
timing <- system.time({
    chr_names <- ArchR:::.availableSeqnames(staging_input, subGroup = "Fragments")
    for (chr in chr_names) {
        f <- ArchR:::.getFragsFromArrow(
            ArrowFile = staging_input, 
            chr = chr, 
            out = "IRanges", 
            cellNames = NULL, 
            method = "fast"
        )
        total_fragments <- total_fragments + length(f)
    }
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    total_fragments = total_fragments
)
write_tsv(results, output_timing)