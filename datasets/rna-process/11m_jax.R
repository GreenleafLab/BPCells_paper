suppressPackageStartupMessages({
    library(BPCells)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_dir <- args[1]
output_dir <- args[2]

h5ad_inputs <- list.files(input_dir) |> sort()

mat_list <- lapply(file.path(input_dir, h5ad_inputs), open_matrix_anndata_hdf5)
mat <- do.call(cbind, mat_list) |>
    convert_matrix_type("uint32_t")
mat <- write_matrix_dir(mat, file.path(output_dir, "bpcells"))

# Need to add feature_metadata for genome to make scanpy happy
write_matrix_10x_hdf5(mat, file.path(output_dir, "10x.h5"), 
    feature_metadata=list("genome"=rep_len("g", nrow(mat))))