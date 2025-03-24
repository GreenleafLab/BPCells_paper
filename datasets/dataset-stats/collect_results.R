suppressPackageStartupMessages({
    library(stringr)
    library(tidyverse)
})

script_dir <- commandArgs(trailingOnly=FALSE)
script_dir <- script_dir[str_detect(script_dir, "--file=")] |> str_replace("--file=", "") |> dirname() |> normalizePath()
script_root <- script_dir
while (!("config_vars.sh" %in% list.files(script_root)) && dirname(script_root) != script_root) {
    script_root <- dirname(script_root)
}

source(file.path(script_root, "config_vars.sh"))
source(file.path(script_dir, "../../utils", "result_collection_utils.R"))

input_path <- file.path(RESULTS_ROOT, "raw/datasets/dataset-stats")
output_path <- file.path(RESULTS_ROOT, "data_tables/datasets/")

dir.create(output_path, recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
rna <- files[str_detect(files, "rna-.*.tsv")]
atac <- files[str_detect(files, "atac-.*.tsv")]

rna <- read_xsv_dataset(rna, read_tsv, "dataset-stats/rna-([^.]+)\\.tsv", c("dataset"),
    show_col_types=FALSE)

atac <- read_xsv_dataset(atac, read_tsv, "dataset-stats/atac-([^.]+)\\.tsv", c("dataset"),
    show_col_types=FALSE)

write_tsv(rna, file.path(output_path, "rna.tsv"))
write_tsv(atac, file.path(output_path, "atac.tsv"))
