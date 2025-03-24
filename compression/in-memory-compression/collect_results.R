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
source(file.path(script_dir, "../bitwidth-helpers.R"))

input_path <- file.path(RESULTS_ROOT, "raw/compression/in-memory-compression")
output_path <- file.path(RESULTS_ROOT, "data_tables/compression/in-memory-compression.tsv.gz")

dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)


files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
bpcells_files <- files[str_detect(files, "bpcells_rep")]
lzbench_files <- files[str_detect(files, "lzbench_rep")]
blosc_files <- files[str_detect(files, "blosc_rep")]


bpcells_raw <- read_xsv_dataset(bpcells_files, read_tsv,
    "in-memory-compression/([^_]+)_(rep[0-9])/([^/]+).tsv",
    c("runner", "replicate_filename", "sample"),
    show_col_types=FALSE
)
stopifnot(all.equal(bpcells_raw$replicate_filename, bpcells_raw$replicate))

bpcells_sizes <- bind_rows(
    matrix_storage(file.path(DATA_ROOT, "rna/130k_thymus_atlas/bpcells")) |> mutate(sample="130k_thymus"),
    matrix_storage(file.path(DATA_ROOT, "rna/130k_thymus_atlas/bpcells_transpose")) |> mutate(sample="130k_thymus_transpose"),
    fragments_storage(file.path(DATA_ROOT, "atac/10k_pbmc/bpcells_filtered/10k_pbmc")) |> mutate(sample="10k_pbmc")

) |> mutate(input_bytes=4*elements) |> select(!c(folder, elements))

bpcells <- bpcells_raw |>
    mutate(
        codec="bpcells",
        field=str_match(expression, "_([a-z]+)")[,2],
        operation=str_match(expression, "(read|write)_")[,2],
        level=NA_integer_,
        software_version=NA,
        filter=NA
    ) |>
    select(!c(expression, `itr/sec`, `gc/sec`, input, replicate_filename, total_time, n_gc, mem_alloc), iterations=n_itr) |>
    pivot_wider(
        names_from="operation", 
        names_glue="{operation}_{.value}",
        values_from=c(min, median, iterations)
    ) |>
    left_join(bpcells_sizes, by=c("field", "sample")) |>
    select(sample, field, codec, level, filter, runner, replicate, write_min, read_min, write_median, read_median, write_iterations, read_iterations, bytes, input_bytes, software_version)

stopifnot(setequal(bpcells$sample, c("10k_pbmc", "130k_thymus", "130k_thymus_transpose")))

# Gather 

lzbench_raw <- read_xsv_dataset(lzbench_files, read_csv,
    "in-memory-compression/([^_]+)_(rep[0-9])/(.+)-(.+)",
    c("runner", "replicate", "field", "sample"),
    show_col_types=FALSE
)
lzbench <- lzbench_raw |>
    mutate(
        codec=str_match(`Compressor name`, "^([^ ]+) ?")[,2],
        software_version=str_match(`Compressor name`, "^[^ ]+ ([0-9.]+)")[,2],
        level=as.integer(str_match(`Compressor name`, "^[^ ]+ [0-9.]+ -(-?[0-9]+)")[,2]),
        write_min=`Compression time in us`/1e6,
        read_min=`Decompression time in us`/1e6,
        bytes=`Compressed size`,
        input_bytes=`Original size`,
        write_median=NA,
        read_median=NA,
        write_iterations=NA,
        read_iterations=NA,
        filter=NA
    ) |>
    select(sample, field, codec, level, filter, runner, replicate, write_min, read_min, write_median, read_median, write_iterations, read_iterations, bytes, input_bytes, software_version)

blosc_raw <- read_xsv_dataset(blosc_files, read_tsv,
    "in-memory-compression/([^_]+)_(rep[0-9])/([^-]+)-([^-]+)-(.+).tsv",
    c("runner", "replicate_filename", "filter_filename", "field", "sample"),
    show_col_types=FALSE
)
stopifnot(all.equal(blosc_raw$replicate, blosc_raw$replicate_filename))
stopifnot(all(blosc_raw$filter == blosc_raw$filter_filename | (blosc_raw$filter=="none" & blosc_raw$method=="numpy")))

blosc <- blosc_raw |>
    select(sample, field, codec, level, filter, runner=method, replicate, write_min, read_min, write_median, read_median, write_iterations, read_iterations, bytes, input_bytes, software_version)

bind_rows(
    bpcells,
    lzbench,
    blosc
) |> write_tsv(output_path)

