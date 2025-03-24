suppressPackageStartupMessages({
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

input_path <- file.path(RESULTS_ROOT, "raw/compression/rna-1M-cell/")
output_path <- file.path(RESULTS_ROOT, "data_tables/compression/rna-1M-cell.tsv")
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)



files <- list.files(input_path, full.names=TRUE, recursive=TRUE) |> normalizePath()
tsv_files <- files[str_detect(files, ".*.tsv")]
gnutime_files <- files[str_detect(files, ".*.gnutime.txt")]

tsv_data <- read_xsv_dataset(tsv_files, read_tsv, 
    "rna-1M-cell/(rep[0-9])/(bpcells|scanpy).tsv",
    c("replicate", "runner"),
    show_col_types=FALSE) |>
    select(!c(runner))

size_lookup <- tsv_data |> distinct(format, size)

gnutime_data <- read_xsv_dataset(gnutime_files, read_gnutime,
    "rna-1M-cell/(rep[0-9])/cp_(10x|h5ad).gnutime.txt",
    c("replicate", "format")) |>
    mutate(tool = "cp", write=seconds, read=NA, read_seconds_cpu=NA, write_seconds_cpu=NA) |>
    inner_join(size_lookup, by=c("format")) |>
    select(!c(seconds, max_rss))



data <- bind_rows(tsv_data, gnutime_data)


write_tsv(data, output_path)

