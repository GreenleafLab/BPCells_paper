library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)

output_dir <- args[1]

#Note: there are 113 cells listed in this as cellclass=GABA, subclass=MBGA, celltype=CTXMIX 
# those cells don't show up anywhere in the bedpe files, so I'll just filter them out here
read_tsv("http://catlas.org/catlas_downloads/humanbrain/Supplementary_Tables/Table%20S3%20%e2%80%93%20Metatable%20and%20annotation%20of%20single%20nuclei.gz") |>
    mutate(cell_id = str_c(sample, ".", barcode)) |>
    filter(celltype != "CTXMIX") |>
    pull(cell_id) |>
    write_lines(file.path(output_dir, "cell-barcodes.txt"))