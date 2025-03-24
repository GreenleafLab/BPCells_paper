suppressPackageStartupMessages({
    library(BPCells)
    library(tidyverse)
    library(presto)
})


args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

tmp_dir <- args[1]
output_timing <- args[2]
output_results <- args[3]


mat <- Seurat::Read10X_h5(file.path(tmp_dir, "input.h5"), use.names=FALSE)
mat <- multiply_cols(mat, 1/colSums(mat))

clusts <- readLines(file.path(tmp_dir, "clusts.txt"))
stopifnot(length(clusts) == ncol(mat))


timing <- system.time({
    res <- wilcoxauc(mat, clusts)
})

write_csv(as_tibble(res), output_results)

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    genes = nrow(mat),
    clusts = length(unique(clusts))
)
write_tsv(results, output_timing)
