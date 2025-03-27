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

input_path <- file.path(RESULTS_ROOT, "raw/cellxgene/03_mean_variance")
if (IS_LAPTOP=="true") {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/03_mean_variance-laptop.tsv")
} else {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/03_mean_variance.tsv")
}
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
timing <- files[str_detect(files, "timing.tsv")]
gnutime <- files[str_detect(files, "timing.gnutime.txt")]


#' Variant to allow additional flexibility with NA-valued groups
#' Read multiple csv files of the same structure,
#' returning a concatenated table with additional columns parsed
#' from the file name
#' @param files list of files
#' @param reader file like read_tsv or read_csv to do the actual reading
#' @param regex regex with groups to parse file names
#' @param group_names column names for the matched groups
#' @param ... Additional arguments to read_csv
read_xsv_dataset_flexible <- function(files, reader, regex, group_names, ...) {
  matches <- str_match(files, regex)
  res <- tibble(file = files)
  for (i in seq_along(group_names)) {
    res[[group_names[i]]] <- matches[,i+1]
    if(any(is.na(matches[,i+1]))) {
      fail_example <- which(is.na(matches[,i+1]))[1]
      warning(sprintf("Got NA group matches in %s for %s lines, e.g. line: %s", group_names[i], sum(is.na(matches[,i+1])), files[fail_example]))
    }
  }
  all_na <- rowMeans(is.na(matches)) == 1
  if (any(all_na)) {
    fail_example <- which(all_na)[1]
    stop(sprintf("Got all NAs for %s lines, e.g. line: %s", sum(all_na), files[fail_example]))
  }
    
  res$data <- lapply(files, reader, ...)
  unnest(select(res, !file), data)
}

runtime <- read_xsv_dataset_flexible(timing, read_tsv,
    "03_mean_variance/(tiledb|bpcells)_(genemajor|cellmajor)?_?(norm|lognorm)?_?([0-9]*)_(rep[0-9]+)/timing\\.tsv",
    c("tool", "format", "normalization", "threads", "replicate"),
    show_col_types=FALSE) 

write_tsv(runtime, output_path)