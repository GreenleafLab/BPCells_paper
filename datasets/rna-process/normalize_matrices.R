suppressPackageStartupMessages({
    library(BPCells)
    library(stringr)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

input_file <- args[1]
output_dir <- args[2]

# Take an input 10x or AnnData hdf5 file and write out to all the required formats
if (str_detect(input_file, "\\.h5$")) {
    input_mat <- open_matrix_10x_hdf5(input_file)
} else if (str_detect(input_file, "\\.h5ad$")) {
    input_mat <- open_matrix_anndata_hdf5(input_file) %>%
        convert_matrix_type("uint32_t")
}

input_mat <- write_matrix_dir(input_mat, file.path(output_dir, "bpcells"))

# Need to add feature_metadata for genome to make scanpy happy
write_matrix_10x_hdf5(input_mat, file.path(output_dir, "10x.h5"), 
    feature_metadata=list("genome"=rep_len("g", nrow(input_mat))))
