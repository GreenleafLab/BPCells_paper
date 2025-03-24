suppressPackageStartupMessages({
    library(BPCells)
    library(Seurat)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 5)

input_dir <- args[1]
tmp_dir <- args[2]
output_var_genes <- args[3]
output_pca <- args[4]
output_clusts <- args[5]

threads <- as.integer(Sys.getenv("OMP_NUM_THREADS"))

cat(sprintf("Variable genes start: %s\n", Sys.time()))
variable_genes <- open_matrix_dir(input_dir) |> 
    CreateSeuratObject() |>
    NormalizeData() |> 
    FindVariableFeatures() |>
    VariableFeatures() 

writeLines(variable_genes, output_var_genes)

cat(sprintf("Normalizing matrix: %s\n", Sys.time()))
# Store the subset sparse normalized matrix to disk
mat <- open_matrix_dir(input_dir)


mat <- mat[rowSums(mat) > 0,] 
mat <- multiply_cols(mat, 10000/colSums(mat))
mat <- log1p(mat[variable_genes,])

gene_stats <- matrix_stats(mat, row_stats="variance", threads=threads)$row_stats
mat <- mat / sqrt(gene_stats["variance",]) 
mat <- write_matrix_dir(mat, file.path(tmp_dir, "normalized_mat"), compress=FALSE)

mat <- mat - gene_stats["mean",] / sqrt(gene_stats["variance",])

cat(sprintf("Running PCA: %s\n", Sys.time()))
gc()
svd <- svds(mat, k=32, threads=threads)

svd[["rownames"]] <- rownames(mat)
svd[["colnames"]] <- colnames(mat)

saveRDS(svd, output_pca, compress=FALSE)


cell_embedding <- t(t(svd$v) * svd$d)
rownames(cell_embedding) <- colnames(mat)

rm(svd)
gc()

cat(sprintf("Running KNN: %s\n", Sys.time()))

knn <- knn_hnsw(cell_embedding, ef=500, threads=threads, k=30)

cat(sprintf("Running Clustering: %s\n", Sys.time()))

# Use geodesic_graph and cluster_leiden to help run faster on large datasets
# Use "modularity" to make cluster count reasonable across dataset sizes
cat(sprintf("Creating graph from knn: %s\n", Sys.time()))
graph <- knn_to_geodesic_graph(knn, threads=threads)
rm(knn)
gc()
cat(sprintf("Runinng leiden clustering: %s\n", Sys.time()))

clusts <- cluster_graph_leiden(graph, resolution=1, objective_function="modularity")
writeLines(as.character(clusts), output_clusts)

cat(sprintf("Done: %s\n", Sys.time()))
