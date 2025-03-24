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

input_path <- file.path(RESULTS_ROOT, "raw/cellxgene/02_matrix_slicing")

if (IS_LAPTOP == "true") {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/02_matrix_slicing-laptop.tsv")
} else {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/02_matrix_slicing.tsv")
}
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
timing <- files[str_detect(files, ".tsv")]
gnutime <- files[str_detect(files, ".gnutime.txt")]

runtime <- read_xsv_dataset(timing, read_tsv,
    "/(tiledb|tiledb_norm|bpcells_[^_]+_?[^_]+)_(rep[0-9]+)/(gene|cell)_([0-9]+)/(random|sequential)_([0-9]+)\\.tsv",
    c("tool_format", "replicate", "axis", "slice_size", "slice_type", "slice_id"),
    show_col_types=FALSE) |>
    mutate(
        tool = str_split_fixed(tool_format, "_", 2)[,1],
        format = str_split_fixed(tool_format, "_", 2)[,2]
    ) |>
    select(tool, format, everything()) |>
    select(-tool_format)

write_tsv(runtime, output_path)