suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(ArchR)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
output_timing <- args[2]
staging_dir <- args[3]
ref_dir <- args[4]

genome <- readLines(file.path(dataset_dir, "genome.txt"))

old_wd <- getwd()
setwd(staging_dir)

addArchRGenome(genome)
proj <- ArchRProject(
    file.path(staging_dir, "sample.arrow"),
    outputDirectory=tempfile(),
    copyArrows=FALSE
)

timing <- system.time({
    addTileMatrix(
        input = proj,
        blacklist = NULL,
        binarize = FALSE,
        excludeChr = c("chrM"),
        threads = 1,
        force=TRUE
    )
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    tool = "archr",
)
saveRDS(proj, file.path(staging_dir, "archr_proj.rds"), compress=FALSE)
write_tsv(results, output_timing)
