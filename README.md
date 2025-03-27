# BPCells Paper code repository

This repository holds the benchmarking code and data tables for the manuscript "Scalable high-performance single cell data analysis with BPCells." The BPCells package itself lives at [github.com/bnprks/BPCells](https://github.com/bnprks/BPCells) and the package documentation is [here](https://bnprks.github.io/BPCells).

## Re-analyzing and plotting data

Data tables from benchmarking live in the [`results/data_tables`](results/data_tables) folder, mostly in tsv format. File paths mirror the structure of the benchmarking folders themselves (see "Figure <-> experiment mapping" below).

R scripts for plotting are in [`results/plots`](results/plots). 

See [`results/README.md`](results/README.md) for details on results file contents and plotting scripts.

## Inspecting benchmark code or re-running benchmarks

Benchmark code lives under the folders `atac-timing`, `cellxgene`, `compression`, `datasets`, and `rna-timing` (see respective README.md files for details). Each individual benchmarking experiment has a single sub-folder. To get started reviewing a particular benchmark, look at the `gen_tasks.py` file and the commands it prints to a `tasks.txt` file (or just skip to the worker scripts based on naming conventions).

Benchmarks are run using an ad hoc system using [`arrayjob/run.py`](arrayjob/run.py), [`config_vars.sh`](config_vars.sh), and per-experiment `gen_tasks.py` files within benchmarking subfolders. See [`arrayjob/README.md`](arrayjob/README.md) for details on re-running benchmarks. Please note that some benchmarks can take a very large amount of compute time to run all replicates and tools. You may want to modify `gen_tasks.py` to reduce the number of datasets, replicates, or tools used for certain benchmarks.

The benchmarking singularity container was converted from a docker image. See [`docker/README.md`](docker/README.md) for details on software versions and where to download the container image.


## Figure <-> experiment mapping

| Figure                 | Experiment path                   |
| ---------------------- | --------------------------------- |
| Fig 1b-d, Fig S1a-c    | rna-timing/pca-benchmark          |
| Fig 1e-g, Fig S1 e+f   | atac-timing/peak-tile-timing      |
| Fig S1d                | rna-timing/marker-genes           |
| Fig 2b-d               | compression/rna-1M-cell           |
| Fig 2e+j, Fig S3a      | compression/in-memory-compression |
| Fig 2g-i               | compression/fragments-read-write  |
| Fig S2                 | compression/bitwidth-stats        |
| Fig S3b                | rna-timing/matrix-transpose       |
| Fig S3c                | atac-timing/merge-fragments       |
| Fig 3a+b, Fig S4a      | cellxgene/01_subset_unique_cells  |
| Fig 3c+d, Fig S4 b,e-g | cellxgene/02_matrix_slicing       |
| Fig 3e, Fig S4 c+d     | cellxgene/03_mean_variance        |
| Fig 3f-h               | cellxgene/04_pca                  |
