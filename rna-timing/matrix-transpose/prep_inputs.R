args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

sample_dir <- args[1]
tmp_dir <- args[2]
tool <- args[3]

stopifnot(tool %in% c("bpcells", "dgCMatrix", "scipy"))

if (tool == "bpcells") {
    system(sprintf("cp -r %s %s", file.path(sample_dir, "bpcells"), file.path(tmp_dir, "input")))
    # Delete the row/col names from BPCells as at high cell counts it negatively affects memory usage
    # and the competing tools aren't being benchmarked on their handling of cell IDs
    system(sprintf("rm %1$s && touch %1$s", file.path(tmp_dir, "input", "row_names")))
    system(sprintf("rm %1$s && touch %1$s", file.path(tmp_dir, "input", "col_names")))
} else {
    file.copy(file.path(sample_dir, "10x.h5"), file.path(tmp_dir, "input.h5"))
}