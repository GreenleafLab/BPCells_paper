suppressPackageStartupMessages({
    library(tidyverse)
})

# Helper functions for calculating compression statistics on BPCells bitpacked files

load_64bit_array <- function(f) {
    x <- readBin(f, "integer", 1e9)
    x <- x[c(-1,-2)]
    hi <- as.numeric(x[c(FALSE, TRUE)])
    lo <- as.numeric(x[c(TRUE, FALSE)])
    hi[hi < 0] <- 2^32 + hi[hi < 0]
    lo[lo < 0] <- 2^32 + lo[lo < 0]
    lo + hi*2^32
}

bitwidth_distribution <- function(folder, field) {
    path <- file.path(folder, paste0(field, "_idx"))
    len <- file.size(path) / 4

    idx <- readBin(path,  "integer", len)[-c(1,2)] |> as.numeric()
    idx[is.na(idx)] <- 2^31 # Handle R interpereting INT_MIN as NA
    idx <- ifelse(idx < 0, idx + 2^32, idx) # We wanted unsigned int, not signed

    # Handle wrap-around for when the 32-bit int overflows
    d <- diff(idx)
    d <- ifelse(d < 0, d + 2^32, d)

    bitwidths <- d * 32/128
    stopifnot(all(bitwidths %in% (0:32)))

    freq <- NULL
    for (i in 0:32) {
        freq <- c(freq, sum(bitwidths == i))
    }
    tibble::tibble(
        folder=folder,
        field=field,
        bitwidth=0:32,
        chunk_count=freq
    )
}

matrix_storage <- function(folder) {
    files <- list.files(folder)
    stopifnot(length(files)==13)
    bytes <- file.size(list.files(folder, full.names=TRUE))
    field <- dplyr::case_when(
        stringr::str_detect(files, "^index_") ~ "index",
        stringr::str_detect(files, "^val_") ~ "val",
        files == "idxptr" ~ "idxptr",
        .default = "metadata"
    )

    # Read 64-bit integers from nonzeros
    nonzeros <- load_64bit_array(file.path(folder, "idxptr"))

    tibble::tibble(
        folder=folder,
        bytes=bytes,
        field=field
    ) %>%
        dplyr::group_by(folder, field) %>%
        dplyr::summarize(bytes=sum(bytes), .groups="drop") %>%
        dplyr::mutate(elements=max(nonzeros))
}

fragments_storage <- function(folder) {
    files <- list.files(folder)
    stopifnot(length(files)==15)
    bytes <- file.size(list.files(folder, full.names=TRUE))
    field <- dplyr::case_when(
        files == "cell_names" ~ "metadata",
        files == "end_max" ~ "end_max",
        stringr::str_detect(files, "^cell_") ~ "cell",
        stringr::str_detect(files, "^start_") ~ "start",
        stringr::str_detect(files, "^end_") ~ "end",
        .default = "metadata"
    )

    nfrags <- load_64bit_array(file.path(folder, "chr_ptr"))

    tibble::tibble(
        folder=folder,
        bytes=bytes,
        field=field
    ) %>%
        dplyr::group_by(folder, field) %>%
        dplyr::summarize(bytes=sum(bytes), .groups="drop") %>%
        dplyr::mutate(elements=max(nfrags))
}

##########################
### HIGH LEVEL HELPERS ###
##########################
# Given a dir and additional metadata via ..., return a tibble
# of data


calculate_bitwidths_fragments <- function(dir) {
    bind_rows(
        bitwidth_distribution(dir, "cell"),
        bitwidth_distribution(dir, "start"),
        bitwidth_distribution(dir, "end")
    ) 
}

# Calculate bitwidths of a matrix given its directory
calculate_bitwidths_matrix <- function(dir) {
    bind_rows(
        bitwidth_distribution(dir, "val"),
        bitwidth_distribution(dir, "index")
    ) 
}

# Get the full set of directories for all matrix variants as a tibble
matrix_paths <- function(original_dir, shuffle_dir) {
    res <- NULL
    for (transpose in c(TRUE, FALSE)) {
        for (ordering in c("original", "shuffle", "mean", "nonzero")) {
            if (!transpose && ordering == "original") {
                path <- file.path(original_dir)
            } else {
                ordering_suffix <- c(
                    "original" = "",
                    "shuffle" = "_shuffle",
                    "mean" = "_sortmean",
                    "nonzero" = "_sortnnz"
                )[ordering]
                transpose_suffix <- ifelse(transpose, "_transpose", "")
                path <- paste0(shuffle_dir, transpose_suffix, ordering_suffix)
            }
            res <- bind_rows(
                res,
                tibble(transpose=transpose, ordering=ordering, folder=path)
            )
        }
    }
    return(res)
}

matrix_value_histogram <- function(mat_dir) {
    h <- open_matrix_dir(mat_dir) |>
        BPCells:::iterate_matrix() |>
        BPCells:::matrix_value_histogram_cpp(1024)
    tibble(
        value = seq_along(h),
        count = h
    )
}