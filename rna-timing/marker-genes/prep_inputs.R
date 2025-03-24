args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)

sample_dir <- args[1]
tmp_dir <- args[2]
tool <- args[3]

stopifnot(tool %in% c("bpcells", "presto", "scanpy", "scanpy_tiecorrect"))

if (tool == "bpcells") {
    system(sprintf("cp -r %s %s", file.path(sample_dir, "bpcells_transpose"), file.path(tmp_dir, "input")))
} else {
    file.copy(file.path(sample_dir, "10x.h5"), file.path(tmp_dir, "input.h5"))
}

file.copy(file.path(sample_dir, "clusts.txt"), file.path(tmp_dir, "clusts.txt"))