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

input_path <- file.path(RESULTS_ROOT, "raw/atac-timing/merge-fragments/")
output_path <- file.path(RESULTS_ROOT, "data_tables/atac-timing/merge-fragments.tsv")
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)



files <- list.files(input_path, full.names=TRUE, recursive=TRUE) |> normalizePath()
tsv_files <- files[str_detect(files, "rep[0-9].tsv")]
gnutime_files <- files[str_detect(files, "rep[0-9].gnutime.txt")]

memory <- read_xsv_dataset(gnutime_files, read_gnutime, 
    "merge-fragments/(unix_sort|bpcells)_([a-z0-9_]+)/(rep[0-9]).gnutime.txt",
    c("tool", "dataset", "replicate"))

time <- read_xsv_dataset(tsv_files, read_tsv,
    "merge-fragments/(unix_sort|bpcells)_([a-z0-9_]+)/(rep[0-9]).tsv", c("tool", "dataset", "replicate"),
    show_col_types=FALSE)

data <- inner_join(time, select(memory, tool, dataset, replicate, max_rss), 
                by=c("tool", "dataset", "replicate"), relationship="one-to-one")

write_tsv(data, output_path)
