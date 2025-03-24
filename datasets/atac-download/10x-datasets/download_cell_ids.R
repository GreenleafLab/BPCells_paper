library(tidyverse)



args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)

output_dir <- file.path(args[1], "atac")


# Mark the cell barcodes
read_csv("https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv") |>
    filter(is__cell_barcode == 1) |>
    pull(barcode) |>
    write_lines(file.path(output_dir, "10k_pbmc/cell-barcodes.txt"))

read_csv("https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_per_barcode_metrics.csv") |>
    filter(is_cell == 1) |>
    pull(barcode) |>
    write_lines(file.path(output_dir, "3k_pbmc/cell-barcodes.txt"))

