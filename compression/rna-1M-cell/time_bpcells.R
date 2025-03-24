suppressPackageStartupMessages({
    library(BPCells)
    library(Matrix)
    library(tidyverse)
})


args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)
path_10x <- args[1]
path_h5ad <- args[2]
scratch_path <- args[3]
output_path <- args[4]

stopifnot(grepl(".h5$", path_10x))
stopifnot(grepl(".h5ad$", path_h5ad))

path_bpcells <- file.path(scratch_path, "bpcells")
path_bpcells_raw <- file.path(scratch_path, "bpcells_mat_raw")
mat1 <- open_matrix_dir(path_bpcells)

# Put this immediately after writing to match the positioning with
# zarr
time_bpcells <- system.time({
    open_matrix_dir(path_bpcells) %>% colSums()
})

# Do write benchmarks starting from in-memory data to avoid simultaneously benchmarking read speed and be consistent with the other tools
mat_mem <- write_matrix_memory(mat1, compress=FALSE)

time_mat_write <- system.time({
    write_matrix_dir(mat_mem, file.path(scratch_path, "bpcells_mat2"))
})

# Do a disk speed control with uncompressed data
time_mat_write_raw <- system.time({
    mat_mem %>% write_matrix_dir(path_bpcells_raw, compress=FALSE)
})
time_bpcells_raw <- system.time({
    open_matrix_dir(path_bpcells_raw) %>% colSums()
})

time_10x <- system.time({
    open_matrix_10x_hdf5(path_10x) %>% colSums()
})

time_h5ad <- system.time({
    open_matrix_anndata_hdf5(path_h5ad) %>% colSums()
})


size_10x <- file.info(path_10x)$size
size_h5ad <- file.info(path_h5ad)$size
size_bpcells <- sum(file.info(list.files(path_bpcells, all.files=TRUE, recursive=TRUE, full.names=TRUE))$size)
size_bpcells_raw <- sum(file.info(list.files(path_bpcells_raw, all.files=TRUE, recursive=TRUE, full.names=TRUE))$size)

data <- tibble(
    tool = "bpcells",
    format = c("10x", "h5ad", "bpcells", "bpcells_raw"),
    read = c(time_10x[3], time_h5ad[3], time_bpcells[3], time_bpcells_raw[3]),
    read_seconds_cpu = c(sum(time_10x[-3]), sum(time_h5ad[-3]), sum(time_bpcells[-3]), sum(time_bpcells_raw[-3])),
    write = c(NA, NA, time_mat_write[3], time_mat_write_raw[3]),
    write_seconds_cpu = c(NA, NA, sum(time_mat_write[-3]), sum(time_mat_write_raw[-3])),
    size = c(size_10x, size_h5ad, size_bpcells, size_bpcells_raw)
)

write_tsv(data, output_path)