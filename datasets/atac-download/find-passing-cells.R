suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 6)

sample_path <- args[1]
output_path <- args[2]
genome <- args[3]
min_frags <- as.numeric(args[4])
min_tss <- as.numeric(args[5])
ref_dir <- args[6]

sample_name <- basename(sample_path)

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))

if (genome == "hg19") {
    genes <- read_gencode_genes(ref_dir, release=19, features="transcript", tags="basic", annotation_set = "comprehensive") %>%
        dplyr::filter(basic)
    blacklist <- read_encode_blacklist(ref_dir, genome="hg19")
} else if (genome == "hg38") {
    genes <- read_gencode_transcripts(ref_dir, release=43, feature="transcript")
    blacklist <- read_encode_blacklist(ref_dir, genome="hg38")
}

qc_stats <- open_fragments_dir(sample_path) |>
    prefix_cell_names(paste0(sample_name, ".")) |>
    select_chromosomes(standard_chrs) |>
    qc_scATAC(genes, blacklist)


qc_stats |>
    dplyr::filter(nFrags >= min_frags, TSSEnrichment >= min_tss) |>
    pull(cellName) |>
    write_lines(output_path)
