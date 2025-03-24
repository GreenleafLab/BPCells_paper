
suppressPackageStartupMessages({
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input <- args[1]
output <- args[2]

readRDS(input) |> 
    rownames() |>
    (\(x) tibble(x=x))() |>
    separate(x, into=c("chr", "start", "end"), sep="_") |>
    write_tsv(output, col_names=FALSE)

