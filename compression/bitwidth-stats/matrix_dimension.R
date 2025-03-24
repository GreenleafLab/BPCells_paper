suppressPackageStartupMessages({
    library(BPCells)
    library(Matrix)
    library(tidyverse)
})

script_dir <- commandArgs(trailingOnly=FALSE)
script_dir <- script_dir[str_detect(script_dir, "--file=")] |> str_replace("--file=", "") |> dirname() |> normalizePath()
source(file.path(script_dir, "../bitwidth-helpers.R"))

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)
input_dir <- args[1]
output_dir <- args[2]

matrix_dims <- readBin(file.path(input_dir, "shape"), integer(), 4)[3:4]

tibble(rows=matrix_dims[1], cols=matrix_dims[2]) |>
    write_tsv(file.path(output_dir, "matrix_shape.tsv"))
