suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(ArchR)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
output_timing <- args[2]
peaks_file <- args[3]
staging_dir <- args[4]

genome <- read_lines(file.path(dataset_dir, "genome.txt"))
old_wd <- getwd()
setwd(staging_dir)

addArchRGenome(genome)
proj <- ArchRProject(
    file.path(staging_dir, "sample.arrow"),
    outputDirectory="archr_project",
    copyArrows=FALSE
)

chrs <- seqnames(getGenomeAnnotation()$chromSizes) |> as.character()

peaks <- read_tsv(peaks_file, col_names=c("chr", "start", "end"), col_types="cii") |>
    mutate(start=start+1) |>
    filter(chr %in% chrs)
stopifnot(nrow(distinct(peaks)) == nrow(peaks))

proj <- addPeakSet(proj, GenomicRanges::makeGRangesFromDataFrame(peaks))

cat(sprintf("Calculating matrices %s\n", Sys.time()))

timing <- system.time({
    addPeakMatrix(
        ArchRProj = proj,
        ceiling = 1e9,
        binarize = FALSE,
        verbose = TRUE,
        threads = 1L,
        parallelParam = NULL,
        force = TRUE
    )
})

cat(sprintf("Outputting timing results %s\n", Sys.time()))

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    tool = "archr",
    peaks = nrow(peaks)
)
write_tsv(results, output_timing)

saveRDS(proj, file.path(staging_dir, "archr_proj.rds"), compress=FALSE)
cat(sprintf("Done %s\n", Sys.time()))