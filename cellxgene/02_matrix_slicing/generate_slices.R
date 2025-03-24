

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

total_count <- as.integer(args[1])
slice_size <- as.integer(args[2])
slice_variants <- as.integer(args[3])
output_dir <- args[4]

# Set seed dependent on slice size and number of variants, since
# otherwise we'll pick the same prefix while sampling
set.seed(11612512L + total_count + slice_size + slice_variants)

# Output random slices
if (slice_size < total_count) {
    for (i in seq_len(slice_variants))  {
        sample.int(n=total_count, size=slice_size) |>
            sort() |>
            as.character() |>
            writeLines(file.path(output_dir, sprintf("random_%03d.txt", i)))
    }
}

# Output sequential slices
for (i in seq_len(slice_variants)) {
    start_coord <- sample.int(n=total_count - slice_size + 1, size=1)
    (seq_len(slice_size) + start_coord - 1) |>
        as.integer() |> 
        as.character() |>
        writeLines(file.path(output_dir, sprintf("sequential_%03d.txt", i)))
}