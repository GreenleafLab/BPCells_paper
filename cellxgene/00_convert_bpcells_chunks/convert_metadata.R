suppressPackageStartupMessages({
    library(tiledbsoma)
    library(tidyverse)
    library(BPCells)
    library(Matrix)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

tiledb_path <- args[1]
ref_dir <- args[2]
output_base <- args[3]

get_table <- function(x) {
    it <- x$read()
    chunks <- list(it$read_next()$to_data_frame())
    while (!it$read_complete()) {
        cat(sprintf("Reading obs: %d, total_rows=%d\n", length(chunks) + 1, sum(as.numeric(lapply(chunks, nrow)))))
        chunks <- c(chunks, list(it$read_next()$to_data_frame()))
    }
    bind_rows(chunks)
}

SOMADataFrameOpen(sprintf("%s/../census_info/datasets", tiledb_path)) %>%
    get_table() %>%
    write_tsv(file.path(output_base, "datasets.tsv"))

### Investigating the cell metadata:
get_obs_table <- function(x) {
    # Loop over rows across all obs (~100-1,000 chunks)
    it <- SOMADataFrameOpen(file.path(tiledb_path, "obs"))$read()
    obs_list <- list(it$read_next()$to_data_frame())
    while (!it$read_complete()) {
        cat(sprintf("Reading obs: %d, total_rows=%d\n", length(obs_list) + 1, sum(as.numeric(lapply(obs_list, nrow)))))
        obs_list <- c(obs_list, list(it$read_next()$to_data_frame()))
    }
    rm(it)
    gc()

    # Check that all the column names match
    stopifnot(all(as.logical(lapply(obs_list, function(x) all.equal(names(x), names(obs_list[[1]]))))))

    all_obs <- list()
    for (n in names(obs_list[[1]])) {
        cat(sprintf("Concatenating column %s\n", n))
        col <- do.call(c, lapply(obs_list, function(df) df[[n]]))
        if (typeof(col) == "character") {
            col <- as.factor(col)
        }
        all_obs[[n]] <- col
        gc()
    }
    all_obs <- as_tibble(all_obs)
    return(all_obs)
}

all_obs <- get_obs_table(x)

# ~4.5 GB file
saveRDS(all_obs, file.path(output_base, "obs.rds"), compress=FALSE)


get_feature_presence <- function(x) {
    f <- SOMASparseNDArrayOpen(file.path(tiledb_path, "ms/RNA/feature_dataset_presence_matrix"))
    it <- f$read()$sparse_matrix()
    feature_presence <- it$read_next()
    count <- 1
    while (!it$read_complete()) {
        cat(sprintf("Reading var: %d, total_nonzeros=%d\n", count, length(feature_presence@i)))
        count <- count + 1
        feature_presence <- feature_presence + it$read_next()
    }
    feature_presence <- as(feature_presence, "nsparseMatrix")
    return(feature_presence)
}

# ~60 MB file, samples x features matrix
feature_presence <- get_feature_presence(x)
saveRDS(feature_presence, file.path(output_base, "feature_presence.rds"), compress=FALSE) 


get_var_list <- function(x) {
    it <- SOMADataFrameOpen(file.path(tiledb_path, "ms/RNA/var"))$read()
    var_list <- list(it$read_next()$to_data_frame())
    while (!it$read_complete()) {
        cat(sprintf("Reading var: %d, total_rows=%d\n", length(var_list) + 1, sum(as.numeric(lapply(var_list, nrow)))))
        var_list <- c(var_list, list(it$read_next()$to_data_frame()))
    }
    rm(it)
    var_table <- bind_rows(var_list)
    return(var_table)
}
var_list <- get_var_list(x)
saveRDS(var_list, file.path(output_base, "var.rds"), compress=FALSE) 



cat("Selecting unique cells\n")
obs_idx <- 1L + all_obs$soma_joinid[all_obs$is_primary_data]
saveRDS(obs_idx, file.path(output_base, "obs_idx.rds"), compress=FALSE)

# List of full gene assays from cellxgene_census_builder: https://github.com/chanzuckerberg/cellxgene-census/blob/bfb40f6acc5b12f31aa24a2f8f244355e876d779/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/globals.py#L328
FULL_GENE_ASSAY = c(
    "EFO:0003755",  # FL-cDNA
    "EFO:0008441",  # full length single cell RNA sequencing
    "EFO:0008747",  # FRISCR
    "EFO:0008763",  # Hi-SCL
    "EFO:0008797",  # MATQ-seq
    "EFO:0008877",  # Quartz-seq
    "EFO:0008930",  # Smart-seq
    "EFO:0008931",  # Smart-seq2
    "EFO:0008956",  # SUPeR-seq
    "EFO:0009810",  # full length single nucleus RNA sequencing
    "EFO:0010004",  # SCRB-seq
    "EFO:0010022",  # Smart-3Seq
    "EFO:0010058",  # Fluidigm C1-based library preparation
    "EFO:0010184",  # Smart-like
    "EFO:0022396",  # TruSeq
    "EFO:0022488",  # Smart-seq3
    "EFO:0022839",  # STORM-seq
    "EFO:0030031",  # SCOPE-chip
    "EFO:0030061",  # mcSCRB-seq
    "EFO:0700016"   # Smart-seq v4
)
stopifnot(all.equal(all_obs$soma_joinid + 1L, seq_len(nrow(all_obs))))
obs_full_gene_mask <- all_obs$assay_ontology_term_id[obs_idx] %in% FULL_GENE_ASSAY
saveRDS(obs_full_gene_mask, file.path(output_base, "obs_full_gene_mask.rds"), compress=FALSE)
cat("Success!\n")
