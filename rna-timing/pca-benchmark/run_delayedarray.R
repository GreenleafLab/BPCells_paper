suppressPackageStartupMessages({
    library(DelayedArray)
    library(DelayedMatrixStats)
    library(HDF5Array)
    library(readr)
    library(tibble)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 5)

sample_dir <- args[1]
data_dir <- args[2]
results_dir <- args[3]
tmp_dir <- args[4]
step <- args[5]

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
if (threads > 1 && step == "normalize") {
    # Can't use parallelism for PCA, because MulticoreParam
    # results in the process hanging for (vec %*% mat) operations
    # (I also tried SnowParam, but it was extremely slow).
    suppressPackageStartupMessages({
        library(BiocParallel)
    })
    setAutoBPPARAM(MulticoreParam(workers=threads))
}

input_matrix <- file.path(sample_dir, "10x.h5")
input_var_genes <- file.path(sample_dir, "variable_genes.txt")

stopifnot(step %in% c("normalize", "pca"))

normalize_path <- file.path(tmp_dir, "normalized.h5")
pca_path <- file.path(data_dir, "pca.rds")

normalize_timing <- file.path(results_dir, "normalize-timing.tsv")
pca_timing <- file.path(results_dir, "pca-timing.tsv")


if (step == "normalize") {
    input_tmp <- file.path(tmp_dir, "input.h5")
    output_tmp <- file.path(tmp_dir, "output.h5")
    file.copy(input_matrix, input_tmp)

    mat <- TENxMatrix(input_tmp)
    variable_genes <- readLines(input_var_genes)
    stopifnot(sum(rownames(mat) %in% variable_genes) == length(variable_genes))
    variable_genes <- rownames(mat)[rownames(mat) %in% variable_genes]

    timing <- system.time({
        # Thanks to Hervé Pagès for providing this edited benchmark code
        # for DelayedArray. Original code commented out below
        setAutoGridMaker("colAutoGrid")
        rs0 <- rowSums(mat)
        cs0 <- colSums(mat)
        mat2 <- mat[rs0 > 0, ]
        mat3 <- t(t(mat2) * 10000 / cs0)
        
        # Simulate as if we were calculating variable genes
        rm3 <- rowMeans(mat3)
        rv3 <- rowVars(mat3, center=rm3)

        mat4 <- log1p(mat3[variable_genes, ])
        mat5 <- mat4 / sqrt(rowVars(mat4))
        writeTENxMatrix(mat5, normalize_path, group="matrix", level=0)

        # Pre-Hervé code version
        # mat <- mat[rowSums(mat) > 0,]
        # mat <- t(t(mat) * 10000 / colSums(mat))
        
        # # Simulate as if we were calculating variable genes
        # invisible(rowMeans(mat))
        # invisible(rowVars(mat))

        # mat <- log1p(mat[variable_genes,])
        # mat <- mat / rowSds(mat)
        # writeTENxMatrix(mat, normalize_path, group="matrix", level=0)
    })
    results <- tibble(
        time_cpu = sum(timing[-3]),
        time_elapsed = timing[3],
        step = step,
        n_ops = NA
    )
    write_tsv(results, normalize_timing)
} else if (step == "pca") {
    mat <- TENxMatrix(normalize_path)
    # We use RSpectra rather than irlba as RSpectra has noticeably lower memory usage
    # and similar time in my experience (it's what BPCells uses internally)
    timing <- system.time({
        row_means <- rowMeans(mat)
        Ax <- function(x, args) {as.numeric(mat %*% x) - row_means * sum(x)}
        Atx <- function(x, args) {as.numeric(x %*% mat) - as.vector(row_means %*% x)}

        svd <- RSpectra::svds(Ax, Atrans=Atx, k=50, dim=dim(mat))
    })

    saveRDS(svd, pca_path, compress=FALSE)

    results <- tibble(
        time_cpu = sum(timing[-3]),
        time_elapsed = timing[3],
        step = step,
        n_ops = svd$nops
    )
    write_tsv(results, pca_timing)
}
