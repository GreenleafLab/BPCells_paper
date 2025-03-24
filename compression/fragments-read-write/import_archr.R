# Import *only* fragments from ArchR, without running QC
suppressPackageStartupMessages({
    library(readr)
    library(tibble)
    library(ArchR)
})


args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4)

dataset_dir <- args[1]
reference_dir <- args[2]
staging_dir <- args[3]
output_timing <- args[4]

genome <- read_lines(file.path(dataset_dir, "genome.txt"))
addArchRGenome(genome)
chrom_sizes <- getArchRGenome(genomeAnnotation=TRUE)$chromSizes

old_wd <- getwd()
setwd(staging_dir)

timing <- system.time({
    tmp <- ArchR:::.tabixToTmp(
        tabixFile = file.path(staging_dir, "sample.fragments.tsv.gz"), 
        tmpFile = file.path(staging_dir, "tmp.arrow"),
        sampleName = "sample", 
        validBC = NULL,
        chromSizes = chrom_sizes, 
        nChunk = 5, # Default to top-level function, though 3 is listed in the lower-level functions
        threads = 1,
        gsubExpression = NULL, 
        prefix = "", 
        verbose = FALSE, 
        tstart = NULL, 
        logFile = NULL
    )

    out <- ArchR:::.tmpToArrow(
        tmpFile = tmp, 
        outArrow = file.path(staging_dir, "sample.arrow"), 
        genome = NULL, # This argument never actually used 
        minFrags = 0, 
        maxFrags = 1e9, 
        sampleName = "sample", 
        prefix = "", 
        threads = 1,
        verbose = FALSE, 
        tstart = NULL, 
        chromSizes = chrom_sizes, 
        removeFilteredCells = FALSE, 
        logFile = NULL
    )
})

results <- tibble(
    time_cpu = sum(timing[-3]),
    time_elapsed = timing[3],
    bytes = file.size(out)
)
write_tsv(results, output_timing)