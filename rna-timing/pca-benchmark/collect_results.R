suppressPackageStartupMessages({
    library(tidyverse)
    library(hdf5r)
})

script_dir <- commandArgs(trailingOnly=FALSE)
script_dir <- script_dir[str_detect(script_dir, "--file=")] |> str_replace("--file=", "") |> dirname() |> normalizePath()
script_root <- script_dir
while (!("config_vars.sh" %in% list.files(script_root)) && dirname(script_root) != script_root) {
    script_root <- dirname(script_root)
}

source(file.path(script_root, "config_vars.sh"))
source(file.path(script_dir, "../../utils", "result_collection_utils.R"))

input_path <- file.path(RESULTS_ROOT, "raw/rna-timing/pca-benchmark")
output_dir <- file.path(RESULTS_ROOT, "data_tables/rna-timing/pca-benchmark")

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
timing <- files[str_detect(files, "-timing.tsv")]
gnutime <- files[str_detect(files, ".gnutime.txt")]
file_size <- files[str_detect(files, "staged-file-size.tsv")]
# Filter out runs that crashed from timeout, etc
gnutime <- gnutime[file.size(gnutime) != 0]

runtime <- read_xsv_dataset(timing, read_tsv,
    "pca-benchmark/([^/]+)/([0-9]+)__([0-9a-z_]+)__(rep[0-9])/(normalize|pca)-timing.tsv",
    c("tool", "threads", "dataset", "replicate", "step_name"),
    show_col_types=FALSE
)

# Had a one-off bug where scanpy/scanpy_dask wrote the wrong step name in their tsv, but the file name is fine
runtime$step <- runtime$step_name

memory <- read_xsv_dataset(gnutime, read_gnutime,
    "pca-benchmark/([^/]+)/([0-9]+)__([0-9a-z_]+)__(rep[0-9])/(normalize|pca).gnutime.txt",
    c("tool", "threads", "dataset", "replicate", "step")
)

# The rule for keeping this is: if PCA timing passed, then the run passed and everything should be present

data <- runtime |>
    inner_join(memory, by=c("tool", "threads", "dataset", "replicate", "step"), relationship="one-to-one") |>
    select(tool, dataset, step, threads, replicate, time_cpu, time_elapsed, max_rss, n_ops) |>
    group_by(tool, dataset, threads, replicate) |>
    filter("pca" %in% step) |>
    ungroup()

write_tsv(data, file.path(output_dir, "performance.tsv"))


file_size <- read_xsv_dataset(file_size, read_tsv,
    "pca-benchmark/([^/]+)/([0-9]+)__([0-9a-z_]+)__(rep[0-9])/staged-file-size.tsv",
    c("tool", "threads", "dataset", "replicate"),
    show_col_types=FALSE)

write_tsv(file_size, file.path(output_dir, "staged_file_size.tsv"))

# Collect PCA comparison results: 
# - All non-randomized comparison at 130k cells
# - Scanpy vs BPCells and scanpy randomized
pca_dir <- file.path(DATA_ROOT, "rna/pca")

read_pcs <- function(method, dataset, replicate) {
    dir <- sprintf("%s/%s/1__%s__%s", pca_dir, method, dataset, replicate)
    file <- list.files(dir, full.names=TRUE)
    if (str_detect(file, ".rds")) {
        pc <- readRDS(file)
        if (is.list(pc)) return(list(cell=pc$v, gene=pc$u))
        else return(list(cell=pc@cell.embeddings, gene=pc@feature.loadings))
    } else {
        h5 <- hdf5r::H5File$new(file)
        return(list(cell=h5[["obsm/X_pca"]][,] |> t(), gene=h5[["varm/PCs"]][,] |> t()))
    }
}

correlate_pcs <- function(ref, alt, dataset, replicate) {
    pc1 <- read_pcs(ref, dataset, replicate)
    pc2 <- read_pcs(alt, dataset, replicate)
    t1 <- cor(pc1$cell, pc2$cell) |>
        cor_to_tibble() |>
        mutate(axis="cell")
    t2 <- cor(pc1$gene, pc2$gene) |>
        cor_to_tibble() |>
        mutate(axis="gene")
    bind_rows(t1, t2)
}

cor_to_tibble <- function(cor_mat) {
    rownames(cor_mat) <- seq_len(nrow(cor_mat))
    colnames(cor_mat) <- seq_len(ncol(cor_mat))
    as_tibble(cor_mat, rownames="PC_ref") |>
        pivot_longer(!PC_ref, names_to="PC_alt", values_to="pearson") |>
        mutate(PC_ref = as.integer(PC_ref), PC_alt=as.integer(PC_alt))
}

non_random_comparisons <- expand_grid(
        ref=c("bpcells_stage_float", "scanpy", "delayedarray", "seurat", "bpcells_stage_int"), 
        alt=c("bpcells_stage_float", "scanpy", "delayedarray", "seurat", "bpcells_stage_int")
    ) |>
    mutate(dataset="130k_thymus_atlas", replicate="rep1") |>
    filter(ref < alt)

random_comparisons <- expand_grid(
        ref=c("scanpy"),
        alt=c("bpcells_randomized", "scanpy_dask", "bpcells_stage_float"),
        replicate=c(paste0("rep", 1:5))
    ) |>
    mutate(dataset="1m_neurons")
comparisons <- bind_rows(non_random_comparisons, random_comparisons) |>
    rowwise() |>
    mutate(data = list(correlate_pcs(ref, alt, dataset, replicate))) |>
    unnest(data)

write_tsv(comparisons, file.path(output_dir, "accuracy.tsv.gz"))