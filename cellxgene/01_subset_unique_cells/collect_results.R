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

input_path <- file.path(RESULTS_ROOT, "raw/cellxgene/01_subset_unique_cells")
if (IS_LAPTOP=="true") {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/01_subset_unique_cells-laptop.tsv")
} else {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/01_subset_unique_cells.tsv")
}
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
timing <- files[str_detect(files, ".tsv")]
gnutime <- files[str_detect(files, ".gnutime.txt")]

runtime <- read_xsv_dataset(timing, read_tsv,
    "/(tiledb|bpcells)_(genemajor|cellmajor|norm|normalized|raw)_?([0-9]*)_(rep[0-9]+)\\.tsv",
    c("tool", "format", "compression_level", "replicate"),
    show_col_types=FALSE)

memory <- read_xsv_dataset(gnutime, read_gnutime,
    "/(tiledb|bpcells)_(genemajor|cellmajor|norm|normalized|raw)_?([0-9]*)_(rep[0-9]+)\\.gnutime\\.txt",
    c("tool", "format", "compression_level", "replicate"))

data <- runtime |>
    inner_join(memory, by=c("tool", "format", "compression_level", "replicate"), relationship="one-to-one") |>
    select(!"seconds")

write_tsv(data, output_path)