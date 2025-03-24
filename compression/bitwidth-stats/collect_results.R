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

input_path <- file.path(RESULTS_ROOT, "raw/compression/bitwidth-stats/")
output_dir <- file.path(RESULTS_ROOT, "data_tables/compression/bitwidth-stats")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)



files <- list.files(input_path, full.names=TRUE, recursive=TRUE) |> normalizePath()

matrix_bitwidth_files <- files[str_detect(files, ".*matrices.*bitwidths.tsv")]
fragment_bitwidth_files <- files[str_detect(files, ".*fragments.*bitwidths.tsv")]

matrix_storage_files <- files[str_detect(files, ".*matrices.*storage.tsv")]
fragment_storage_files <- files[str_detect(files, ".*fragments.*storage.tsv")]

value_histogram_files <- files[str_detect(files, "value_histogram.tsv")]
matrix_shape_files <- files[str_detect(files, "matrix_shape.tsv")]

matrix_bitwidth_raw <- read_xsv_dataset(matrix_bitwidth_files, read_tsv, 
    "matrices/([^/]+)__([^/]+)/bitwidths.tsv",
    c("dataset", "matrix"),
    show_col_types=FALSE)

stopifnot(setequal(c("rna", "rna_transpose", "peaks", "peaks_transpose", "tiles", "tiles_transpose"), matrix_bitwidth_raw$matrix))
matrix_bitwidth <- matrix_bitwidth_raw |>
    mutate(
        cell_major = matrix %in% c("rna", "peaks_transpose", "tiles_transpose"),
        matrix = str_remove(matrix, "_transpose")
    ) |>
    select(dataset, matrix, cell_major, ordering, everything())
write_tsv(matrix_bitwidth, file.path(output_dir, "matrix_bitwidth.tsv.gz"))

fragment_bitwidth_raw <- read_xsv_dataset(fragment_bitwidth_files, read_tsv, 
    "fragments/([^/]+)/bitwidths.tsv",
    c("dataset"),
    show_col_types=FALSE)

fragment_bitwidth <- fragment_bitwidth_raw |>
    mutate(
        filtered = str_detect(dataset, "_filtered"),
        dataset = str_remove(dataset, "_filtered")
    ) |>
    select(dataset, sample, filtered, everything())
write_tsv(fragment_bitwidth, file.path(output_dir, "fragment_bitwidth.tsv.gz"))

matrix_storage_raw <- read_xsv_dataset(matrix_storage_files, read_tsv, 
    "matrices/([^/]+)__([^/]+)/storage.tsv",
    c("dataset", "matrix"),
    show_col_types=FALSE) 


stopifnot(setequal(c("rna", "rna_transpose", "peaks", "peaks_transpose", "tiles", "tiles_transpose"), matrix_storage_raw$matrix))
matrix_storage <- matrix_storage_raw |>
    mutate(
        cell_major = matrix %in% c("rna", "peaks_transpose", "tiles_transpose"),
        matrix = str_remove(matrix, "_transpose")
    ) |>
    select(dataset, matrix, cell_major, ordering, everything())
write_tsv(matrix_storage, file.path(output_dir, "matrix_storage.tsv"))



fragment_storage_raw <- read_xsv_dataset(fragment_storage_files, read_tsv, 
    "fragments/([^/]+)/storage.tsv",
    c("dataset"),
    show_col_types=FALSE)

fragment_storage <- fragment_storage_raw |>
    mutate(
        filtered = str_detect(dataset, "_filtered"),
        dataset = str_remove(dataset, "_filtered")
    ) |>
    select(dataset, sample, filtered, everything())
write_tsv(fragment_storage, file.path(output_dir, "fragment_storage.tsv"))



matrix_value_histogram_raw <- read_xsv_dataset(value_histogram_files, read_tsv, 
    "matrices/([^/]+)__([^/]+)/value_histogram.tsv",
    c("dataset", "matrix"),
    show_col_types=FALSE)

stopifnot(setequal(c("rna", "rna_transpose", "peaks", "peaks_transpose", "tiles", "tiles_transpose"), matrix_value_histogram_raw$matrix))
matrix_value_histogram <- matrix_value_histogram_raw |>
    mutate(
        cell_major = matrix %in% c("rna", "peaks_transpose", "tiles_transpose"),
        matrix = str_remove(matrix, "_transpose")
    ) |>
    filter(cell_major) |> select(!cell_major) |>
    select(dataset, matrix, everything())
write_tsv(matrix_value_histogram, file.path(output_dir, "matrix_value_histogram.tsv.gz"))


matrix_shape_raw <- read_xsv_dataset(matrix_shape_files, read_tsv, 
    "matrices/([^/]+)__([^/]+)/matrix_shape.tsv",
    c("dataset", "matrix"),
    show_col_types=FALSE)

stopifnot(setequal(c("rna", "rna_transpose", "peaks", "peaks_transpose", "tiles", "tiles_transpose"), matrix_shape_raw$matrix))
matrix_shape <- matrix_shape_raw |>
    mutate(
        cell_major = matrix %in% c("rna", "peaks_transpose", "tiles_transpose"),
        matrix = str_remove(matrix, "_transpose")
    ) |>
    filter(cell_major) |> select(!cell_major) |>
    select(dataset, matrix, everything())
write_tsv(matrix_shape, file.path(output_dir, "matrix_shape.tsv"))


