suppressPackageStartupMessages({
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

input_path <- file.path(RESULTS_ROOT, "raw/rna-timing/marker-genes/")
output_path <- file.path(RESULTS_ROOT, "data_tables/rna-timing/marker-genes")
dir.create(output_path, recursive=TRUE, showWarnings=FALSE)



files <- list.files(input_path, full.names=TRUE, recursive=TRUE) |> normalizePath()
tsv_files <- files[str_detect(files, "rep[0-9].tsv")]
# Emtpy files from crashed jobs casue issues, so filter those out
gnutime_files <- files[str_detect(files, "rep[0-9].gnutime.txt") & file.size(files) > 0]

memory <- read_xsv_dataset(gnutime_files, read_gnutime, 
    "marker-genes/(bpcells|presto|scanpy|scanpy_tiecorrect)_([0-9][a-z0-9_]+)/(rep[0-9]).gnutime.txt",
    c("tool", "dataset", "replicate"))

time <- read_xsv_dataset(tsv_files, read_tsv,
        "marker-genes/(bpcells|presto|scanpy|scanpy_tiecorrect)_([0-9][a-z0-9_]+)/(rep[0-9]).tsv", 
        c("tool", "dataset", "replicate"),
    show_col_types=FALSE)

data <- inner_join(time, select(memory, tool, dataset, replicate, max_rss), 
                by=c("tool", "dataset", "replicate"), relationship="one-to-one")


# Check that all the finished
replicate_count <- data |>
    group_by(tool, dataset) |>
    summarize(count=n(), .groups="drop") |>
    pull(count)
stopifnot(length(unique(replicate_count)) == 1)

write_tsv(data, file.path(output_path, "performance.tsv"))


# Measure and check consistency of results
cat("Calculating accuracy and correlation of results\n")

read_result_bpcells <- function(path) {
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types=FALSE, progress=FALSE) |>
        select(clust=foreground, gene=feature, pval=p_val_raw) |>
        arrange(clust, gene)
}
read_result_presto <- function(path) {
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types=FALSE, progress=FALSE) |>
        select(clust=group, gene=feature, pval=pval) |>
        arrange(clust, gene)
}
read_result_scanpy <- function(path) {
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types=FALSE, progress=FALSE) |>
        select(clust=clust, gene=names, pval=pvals) |>
        arrange(clust, gene)
}

result_path <- function(tool, dataset, replicate) {
    file.path(DATA_ROOT, sprintf("rna/marker-genes/%s/%s_%s.csv", tool, dataset, replicate))
}

all_comparisons <- NULL
for (dataset in unique(data$dataset)) {
    cat(sprintf("Checking accuracy of %s\n", dataset))
    bpcells_data <- result_path("bpcells", dataset, "rep1") |>
        read_result_bpcells()
    presto_data <- result_path("presto", dataset, "rep1") |>
        read_result_presto()

    scanpy_tiecorrect <- result_path("scanpy_tiecorrect", dataset, "rep1") |>
        read_result_scanpy()
    scanpy <- result_path("scanpy", dataset, "rep1") |>
        read_result_scanpy()

    # Confirm that the other replicates are identical
    for (rep in c("rep2", "rep3", "rep4", "rep5")) {
        r2 <- result_path("bpcells", dataset, "rep2") |> read_result_bpcells()
        if (!is.null(r2)) stopifnot(all.equal(bpcells_data, r2))

        r2 <- result_path("presto", dataset, "rep2") |> read_result_presto()
        if (!is.null(r2)) stopifnot(all.equal(presto_data, r2))

        r2 <- result_path("scanpy_tiecorrect", dataset, "rep1") |> read_result_scanpy()
        if (!is.null(r2)) stopifnot(all.equal(scanpy_tiecorrect, r2))

        r2 <- result_path("scanpy", dataset, "rep1") |> read_result_scanpy()
        if (!is.null(r2)) stopifnot(all.equal(scanpy, r2))
    }

    results <- list(
        bpcells=bpcells_data, 
        presto=presto_data,
        scanpy_tiecorrect=scanpy_tiecorrect,
        scanpy=scanpy
    )

    # Skip comparisons if bpcells was the only tool to complete successfully
    if (3 == sum(as.logical(lapply(results, is.null)))) {
        next
    }

    comparisons <- expand_grid(
        tool1 = names(results),
        tool2 = names(results)
    ) |>
        filter(tool1 < tool2) |>
        rowwise() |>
        filter(!is.null(results[[tool1]]), !is.null(results[[tool2]])) |>
        mutate(
            cor = cor(
                pmax(-290, log10(results[[tool1]]$pval)),
                pmax(-290, log10(results[[tool2]]$pval))
            ),
            mean_abs_diff = mean(abs(
                pmax(-290, log10(results[[tool1]]$pval)) -
                pmax(-290, log10(results[[tool2]]$pval))
            )),
            dataset=dataset
        ) |>
        ungroup()
    all_comparisons <- bind_rows(all_comparisons, comparisons)
}
write_tsv(all_comparisons, file.path(output_path, "accuracy.tsv"))

