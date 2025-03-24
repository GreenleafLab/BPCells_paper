# Convert fragment file to ArchR

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(ArchR)
})


get_abspath <- function(path) {
    file.path(normalizePath(dirname(path)), basename(path))
}

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

dataset_dir <- get_abspath(args[1])
sample_id <- args[2]

input_file <- file.path(dataset_dir, "fragments", paste0(sample_id, ".fragments.tsv.gz"))
output_file <- file.path(dataset_dir, "archr", paste0(sample_id, ".arrow"))
genome <- readLines(file.path(dataset_dir, "genome.txt"))

addArchRGenome(genome)

# Make sure these come out as absolute paths since we're going to change directories
# during actual Arrow creation
input_file <- normalizePath(input_file)

# Change to an actual tempdir so that we can have a more reliable temporary filesystem on OAK
original_dir <- getwd()
tmp <- tempdir()
setwd(tmp)

arrow_res <- createArrowFiles(
    inputFiles = input_file,
    sampleName = sample_id,
    outputNames = str_remove(output_file, fixed(".arrow")),
    minTSS = 0,
    minFrags = 0,
    maxFrags = Inf,
    addTileMat = FALSE,
    addGeneScoreMat = FALSE,
    threads = 1
)

stopifnot(!is.null(arrow_res))
# Filter the arrow file to just the keeper cells

output_filtered <- file.path(dataset_dir, "archr", paste0(sample_id, ".filtered.arrow"))

prefix_names <- "true" == read_lines(file.path(dataset_dir, "cell-barcodes-prefix.txt"))
keeper_cells <- readLines(file.path(dataset_dir, "cell-barcodes.txt"))
arrow_cells <- ArchR:::.availableCells(arrow_res)

if (prefix_names) {
    keeper_cells <- str_c(sample_id, "#", str_remove(keeper_cells, fixed(str_c(sample_id, "."))))
} else {
    keeper_cells <- str_c(sample_id, "#", keeper_cells)
}

keeper_cells <- intersect(keeper_cells, arrow_cells)

ArchR:::.copyArrowSingle(arrow_res, output_filtered, keeper_cells)

stopifnot(
    length(intersect(ArchR:::.availableCells(output_filtered), keeper_cells)) 
    == length(keeper_cells)
)

setwd(original_dir)



