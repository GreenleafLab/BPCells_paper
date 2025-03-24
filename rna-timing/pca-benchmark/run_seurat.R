suppressPackageStartupMessages({
    library(Seurat)
    library(tibble)
    library(readr)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 5)

sample_dir <- args[1]
data_dir <- args[2]
results_dir <- args[3]
tmp_dir <- args[4]
step <- args[5]

input_matrix <- file.path(sample_dir, "10x.h5")
input_var_genes <- file.path(sample_dir, "variable_genes.txt")

stopifnot(step %in% c("normalize", "pca"))

normalized_rds <- file.path(tmp_dir, "normalize-mat.rds")
pca_rds <- file.path(data_dir, "pca.rds")

normalize_timing <- file.path(results_dir, "normalize-timing.tsv")
pca_timing <- file.path(results_dir, "pca-timing.tsv")

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
if (threads > 1) {
    # This seems to only effect normalization, but we include it for completeness
    library(future)
    options(future.globals.maxSize=256e09)
    plan("multicore", workers=threads)
}

if (step == "normalize") {
    s <- Read10X_h5(input_matrix, use.names=FALSE)
 
    variable_genes <- readLines(input_var_genes)
    stopifnot(sum(rownames(s) %in% variable_genes) == length(variable_genes))
    variable_genes <- rownames(s)[rownames(s) %in% variable_genes]


    timing <- system.time({
        s <- CreateSeuratObject(s, min.cells=1)
        s <- NormalizeData(s)

        # Dummy calculation of gene mean + variance
        means <- rowMeans(s)
        invisible(Seurat:::SparseRowVar2(s@assays$RNA@layers$data, means, display_progress = FALSE))
        rm(means)
        VariableFeatures(s) <- variable_genes

        s <- ScaleData(s, scale.max=Inf)
    })
    saveRDS(s, normalized_rds, compress=FALSE)

    results <- tibble(
        time_cpu = sum(timing[-3]),
        time_elapsed = timing[3],
        step = step,
        n_ops = NA
    )
    write_tsv(results, normalize_timing)
} else if (step == "pca") {
    s <- readRDS(normalized_rds)

    timing <- system.time({
        s <- RunPCA(s)
    })

    saveRDS(s@reductions$pca, pca_rds, compress=FALSE)

    results <- tibble(
        time_cpu = sum(timing[-3]),
        time_elapsed = timing[3],
        step = step,
        n_ops = NA
    )
    write_tsv(results, pca_timing)
}
