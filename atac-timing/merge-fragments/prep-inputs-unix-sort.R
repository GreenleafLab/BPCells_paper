suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

project_dir <- args[1]
staging_dir <- args[2]
threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))

dataset <- basename(project_dir)
fragment_files <- Sys.glob(file.path(project_dir, "fragments", "*.fragments.tsv.gz"))
samples <- str_remove(basename(fragment_files), ".fragments.tsv.gz")
fragment_dirs <- file.path(project_dir, "bpcells_filtered", samples)

parallel::mclapply(mc.cores=threads, mc.preschedule=FALSE, seq_along(samples), function(i) {
    open_fragments_dir(fragment_dirs[i]) |>
        write_fragments_10x(file.path(staging_dir, sprintf("%s.fragments.tsv", samples[i])))
})
