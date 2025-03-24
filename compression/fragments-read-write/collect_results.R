suppressPackageStartupMessages({
    library(stringr)
})

script_dir <- commandArgs(trailingOnly=FALSE)
script_dir <- script_dir[str_detect(script_dir, "--file=")] |> str_replace("--file=", "") |> dirname() |> normalizePath()
script_root <- script_dir
while (!("config_vars.sh" %in% list.files(script_root)) && dirname(script_root) != script_root) {
    script_root <- dirname(script_root)
}

source(file.path(script_root, "config_vars.sh"))
source(file.path(script_dir, "../../utils", "result_collection_utils.R"))

input_path <- file.path(RESULTS_ROOT, "raw/compression/fragments-read-write/")
output_dir <- file.path(RESULTS_ROOT, "data_tables/compression/fragments-read-write/")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE) |> normalizePath()

import_timing <- files[str_detect(files, "import.tsv")]
read_timing <- files[str_detect(files, "read.tsv")]
gnutime <- files[str_detect(files, ".gnutime.txt")]

gnutime_table <- read_xsv_dataset(gnutime, read_gnutime,
    "fragments-read-write/([^/]+)/([^.]+)\\.([a-z0-9]+)\\.(rep[0-9]+)\\.(read|import)\\.gnutime\\.txt",
    c("dataset", "sample", "tool", "replicate", "step"))

gnutime_table <- gnutime_table |> rename(gnutime_elapsed=seconds)

import_table <- read_xsv_dataset(import_timing, read_tsv,
    "fragments-read-write/([^/]+)/([^.]+)\\.([a-z0-9]+)\\.(rep[0-9]+)\\.import.tsv",
    c("dataset",  "sample", "tool", "replicate"),
    progress=FALSE,
    show_col_types=FALSE)

read_table <- read_xsv_dataset(read_timing, read_tsv,
    "fragments-read-write/([^/]+)/([^.]+)\\.([a-z0-9]+)\\.(rep[0-9]+)\\.read.tsv",
    c("dataset",  "sample", "tool", "replicate"),
    progress=FALSE,
    show_col_types=FALSE)


gnutime_table |>
    filter(step == "read") |> select(-c(step)) |>
    inner_join(read_table, by=c("dataset", "sample", "tool", "replicate"), relationship="one-to-one") |>
    write_tsv(file.path(output_dir, "read.tsv.gz"))

gnutime_table |>
    filter(step == "import") |> select(-c(step)) |>
    inner_join(import_table, by=c("dataset", "sample", "tool", "replicate"), relationship="one-to-one") |>
    write_tsv(file.path(output_dir, "import.tsv.gz"))


