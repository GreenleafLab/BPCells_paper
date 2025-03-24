suppressPackageStartupMessages({
    library(BPCells)
    library(dplyr)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

input_dir <- args[1]
output_dir <- args[2]
tmp_dir <- args[3]

rds_inputs <- list.files(input_dir) |> sort()
rds_inputs <- rds_inputs[rds_inputs != "df_cell.RDS"] # Not a matrix
rds_inputs <- rds_inputs[rds_inputs != "293T_3T3_gene_count.RDS"] # No cells to retain in final dataset

valid_cells <- file.path(input_dir, "df_cell.RDS") |>
    readRDS() |>
    filter(!is.na(Main_cluster_name)) |>
    pull(sample)

for (rds in rds_inputs) {
    cat(sprintf("Reading %s: %s\n", rds, Sys.time()))
    mat <- readRDS(file.path(input_dir, rds)) |>
        as("IterableMatrix") |>
        convert_matrix_type("uint32_t")
    mat_subset <- mat[,colnames(mat) %in% valid_cells]
    write_matrix_dir(mat_subset, file.path(tmp_dir, rds))
    rm(mat, mat_subset)
    gc()
}

cat(sprintf("Finished reading inputs: %s\n", Sys.time()))

mat_list <- lapply(file.path(tmp_dir, rds_inputs), open_matrix_dir)
mat <- do.call(cbind, mat_list)
mat <- write_matrix_dir(mat, file.path(output_dir, "bpcells"))

# Need to add feature_metadata for genome to make scanpy happy
write_matrix_10x_hdf5(mat, file.path(output_dir, "10x.h5"), 
    feature_metadata=list("genome"=rep_len("g", nrow(mat))))


lapply(file.path(tmp_dir, rds_inputs), unlink, recursive=TRUE)