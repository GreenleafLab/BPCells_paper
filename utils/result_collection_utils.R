
# Helper functions for `collect_results.R` files

suppressPackageStartupMessages({
    library(tidyverse)
})

read_gnutime <- function(file) {
    max_rss <- file %>%
        readChar(2000) %>%
        str_match("Maximum resident set size \\(kbytes\\): (.+)") %>%
        {as.numeric(.[,2])*1024}
    elapsed_time_vec <- file %>%
        readChar(2000) %>%
        str_match("Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): ([0-9.:]+)\n") %>%
        .[,2] %>%
        str_split_1(":") %>%
        as.numeric()
    stopifnot(length(elapsed_time_vec) %in% c(2,3))
    seconds <- 0
    for (t in elapsed_time_vec) {
    seconds <- 60 * seconds + t
    }
    {tibble(max_rss = max_rss, seconds=seconds)}
}

#' Read multiple csv files of the same structure,
#' returning a concatenated table with additional columns parsed
#' from the file name
#' @param files list of files
#' @param reader file like read_tsv or read_csv to do the actual reading
#' @param regex regex with groups to parse file names
#' @param group_names column names for the matched groups
#' @param ... Additional arguments to read_csv
read_xsv_dataset <- function(files, reader, regex, group_names, ...) {
  matches <- str_match(files, regex)
  res <- tibble(file = files)
  for (i in seq_along(group_names)) {
    res[[group_names[i]]] <- matches[,i+1]
    if(any(is.na(matches[,i+1]))) {
      fail_example <- which(is.na(matches[,i+1]))[1]
      stop(sprintf("Failed match regex for %s lines, e.g. line: %s", sum(is.na(matches[,i+1])), files[fail_example]))
    }
  }
  res$data <- lapply(files, reader, ...)
  unnest(select(res, !file), data)
}