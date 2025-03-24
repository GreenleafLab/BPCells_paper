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

input_path <- file.path(RESULTS_ROOT, "raw/atac-timing/peak-tile-timing/")
output_path <- file.path(RESULTS_ROOT, "data_tables/atac-timing/peak-tile-timing.tsv.gz")
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
timing <- files[str_detect(files, ".tsv")]
totals <- files[str_detect(files, ".matrix_sum.txt")]
gnutime <- files[str_detect(files, ".gnutime.txt")]

memory <- read_xsv_dataset(gnutime, read_gnutime, 
    "/([^/]+)/(peaks-10*|tiles)/([^.]+)\\.([a-z0-9]+)\\.(rep[0-9]+)\\.gnutime\\.txt",
    c("dataset", "region_type", "sample", "tool", "replicate"))

runtime <- read_xsv_dataset(timing, read_tsv,
    "/([^/]+)/(peaks-10*|tiles)/([^.]+)\\.([a-z0-9]+)\\.(rep[0-9]+)\\.tsv",
    c("dataset", "region_type", "sample", "tool_filename", "replicate"),
    progress=FALSE,
    show_col_types=FALSE)

stopifnot(all.equal(runtime$tool, runtime$tool_filename))
runtime <- select(runtime, -c("tool_filename"))

matrix_sum <- read_xsv_dataset(totals, function(x) {tibble(matrix_sum=as.numeric(read_lines(x)))},
    "/([^/]+)/(peaks-10*|tiles)/([^.]+)\\.([a-z0-9]+)\\.(rep[0-9]+)\\.matrix_sum\\.txt",
    c("dataset", "region_type", "sample", "tool", "replicate"))


data <- runtime |>
    inner_join(matrix_sum, by=c("dataset", "region_type", "sample", "tool", "replicate"), relationship="one-to-one") |>
    inner_join(rename(memory, gnutime_elapsed=seconds), 
                by=c("dataset", "region_type", "sample", "tool", "replicate"),
                relationship="one-to-one") |>
    select(dataset, sample, tool, region_type, replicate, max_rss, time_cpu, time_elapsed, gnutime_elapsed, matrix_sum, peak_count=peaks)

write_tsv(data, output_path)
