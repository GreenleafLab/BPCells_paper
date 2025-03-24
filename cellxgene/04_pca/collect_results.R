suppressPackageStartupMessages({
    library(stringr)
    library(tidyverse)
})

script_dir <- commandArgs(trailingOnly=FALSE)
script_dir <- script_dir[str_detect(script_dir, "--file=")] |> str_replace("--file=", "") |> dirname() |> normalizePath()
script_root <- script_dir
while (!("config_vars.sh" %in% list.files(script_root)) && dirname(script_root) != script_root) {
    script_root <- dirname(script_root)
}

source(file.path(script_root, "config_vars.sh"))
source(file.path(script_dir, "../../utils", "result_collection_utils.R"))

input_path <- file.path(RESULTS_ROOT, "raw/cellxgene/04_pca")
if (IS_LAPTOP=="true") {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/04_pca-laptop")
} else {
    output_path <- file.path(RESULTS_ROOT, "data_tables/cellxgene/04_pca")
}
dir.create(output_path, recursive=TRUE, showWarnings=FALSE)

files <- list.files(input_path, full.names=TRUE, recursive=TRUE)
timing <- files[str_detect(files, ".tsv")]
gnutime <- files[str_detect(files, ".gnutime.txt")]

runtime <- read_xsv_dataset(timing, read_tsv, 
    "(compress|uncompress)_(exact|randomized)_(rep[0-9]+)\\.tsv",
    c("input", "method", "replicate"),
    show_col_types=FALSE)

memory <- read_xsv_dataset(gnutime, read_gnutime,
    "(compress|uncompress)_(exact|randomized)_(rep[0-9]+)\\.gnutime.txt",
    c("input", "method", "replicate"))

data <- runtime |>
    inner_join(memory, by=c("input", "method", "replicate"), relationship="one-to-one") |>
    select(input, method, replicate, time_cpu, time_elapsed, max_rss)

write_tsv(data, file.path(output_path, "timing.tsv"))

cat(sprintf("Measuring similarity between PC components\n"))

pca_dir <- file.path(DATA_ROOT, "cellxgene-census/2024-07-01/pca")

x <- readRDS(file.path(pca_dir, "compress_exact_rep1.rds"))

cor_to_tibble <- function(cor_mat) {
    rownames(cor_mat) <- seq_len(nrow(cor_mat))
    colnames(cor_mat) <- seq_len(ncol(cor_mat))
    as_tibble(cor_mat, rownames="PC_ref") |>
        pivot_longer(!PC_ref, names_to="PC_alt", values_to="pearson") |>
        mutate(PC_ref = as.integer(PC_ref), PC_alt=as.integer(PC_alt))
}

cor_tables <- list()
for (alt in c("compress_randomized_rep1", "uncompress_exact_rep1", "uncompress_randomized_rep1")) {
    cat(sprintf("Running correlation with %s\n", alt))
    y <- readRDS(file.path(pca_dir, paste0(alt, ".rds")))
    gc()
    t1 <- cor(x$u, y$u) |> 
        cor_to_tibble() |>
        mutate(
            ref="compress_exact_rep1",
            alt=alt,
            axis="gene",
            axis_len=nrow(x$u)
        )
    t2 <- cor(x$v, y$v) |> 
        cor_to_tibble() |>
        mutate(
            ref="compress_exact_rep1",
            alt=alt,
            axis="cell",
            axis_len=nrow(x$v)
        )
    cor_tables[[alt]] <- list(bind_rows(t1, t2))
}

write_tsv(bind_rows(cor_tables), file.path(output_path, "pca-cor.tsv.gz"))