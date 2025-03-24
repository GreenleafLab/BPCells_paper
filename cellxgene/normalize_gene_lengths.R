

# For cells that come from full-length assays, divide each row by the feature length
#
# Note: in newer BPCells versions (post 0.3.0), this should be possible to replace with a single line:
# mat[,full_gene_mask] <- multiply_rows(mat[,full_gene_mask], 1/feature_length)
#
# When locking in v0.3.0, however, this needs to be done more manually:
normalize_gene_lengths <- function(matrix, full_gene_mask, feature_length) {
    matrix_list <- if(.hasSlot(matrix, "matrix_list")) matrix@matrix_list else list(matrix)
    chunk_size <- ncol(matrix_list[[1]])

    # Find range of column indices which have full gene assays
    runs <- rle(full_gene_mask)
    range_start <- 1 + cumsum(c(0, runs$lengths))[-length(runs$length)-1][runs$values]
    range_end <- cumsum(c(0, runs$lengths))[-1][runs$values]


    new_mat_list <- list()
    range_idx <- 1L
    for (j in seq_along(matrix_list)) {
        mat <- matrix_list[[j]]
        stopifnot(j == length(matrix_list) || ncol(mat) == chunk_size)
        # Split off ranges that overlap full gene assays
        start_idx <- 1 + chunk_size*(j-1)
        end_idx <- start_idx + ncol(mat) - 1
        if (range_idx > length(range_start) || range_start[range_idx] > end_idx) {
            new_mat_list <- c(new_mat_list, mat)
            next
        }
        
        has_remainder <- TRUE
        while (range_idx <= length(range_start) && range_start[range_idx] <= end_idx && has_remainder) {
            # Take chunk up until the range start
            if (range_start[range_idx] > start_idx) {
                m <- mat[,1:(range_start[range_idx] - start_idx)]
                new_mat_list <- c(new_mat_list, m)
                mat <- mat[,(range_start[range_idx]-start_idx+1):ncol(mat)]
                start_idx <- range_start[range_idx]
            }
            # Take as much of the range as possible
            m <- mat[,1:min(ncol(mat), range_end[range_idx] - start_idx + 1)]
            m <- multiply_rows(m, 1/feature_length)
            new_mat_list <- c(new_mat_list, m)
            has_remainder <- range_end[range_idx] < end_idx
            if (has_remainder) {
                # Resize the remainder before searching for the next range
                mat <- mat[,(range_end[range_idx]-start_idx+2):ncol(mat)]
                start_idx <- range_end[range_idx]+1
            }
            if (range_end[range_idx] <= end_idx) {
                range_idx <- range_idx + 1
            }
        }
        if (has_remainder) {
            new_mat_list <- c(new_mat_list, mat)
        }
        
    }
    if(.hasSlot(matrix, "matrix_list")) {
        matrix@matrix_list <- new_mat_list
    } else {
        matrix <- do.call(cbind, new_mat_list)
    }
    return(matrix)
}
