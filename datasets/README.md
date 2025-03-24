# Dataset prep code

This folder contains code to download RNA + ATAC datasets and create uniform input files for later benchmarking. (Excludes CELLxGENE download). It doesn't produce any specific figures, but these jobs collect the main data dependencies for other benchmarks. 

Expect at least ~1.5TB of disk space to store all the RNA + ATAC data once it is downloaded and processed into copies with various file formats.

- `datasets/atac-download`: 
    - Download ATAC datasets, create BPCells unfiltered fragment files, and call cell barcodes passing quality filters
    - Produces:
        - `{DATA_ROOT}/atac/[dataset]/` 
            - `fragments/[sample].fragments.tsv.gz`: 10x-format fragment files
            - `bpcells/[sample]/`: BPCells-format fragment files
            - `cell-barcodes.txt`, `cell-barcodes/[sample].txt`: Cell barcodes passing quality filtering, or a directory of per-sample barcodes for multi-sample datasets (which will be consolidated into a single list during `atac-process`).
            - `peaks.bed`: Bed-format peak set
            - `genome.txt`: Genome used for the dataset
            - `cell-barcodes-prefix.txt`: "true" or "false" if cell barcodes need to be prefixed by sample names
        - `{DATA_ROOT}/reference_metadata`: Genome reference annotation files
    - Task summary: 209/209 completed. Total time: 1 day, 23:03:07
- `datasets/atac-process`: 
    - Process downloaded ATAC files into uniform inputs
    - Requires: `atac-download`
    - Produces:
        - `{DATA_ROOT}/atac/[dataset]/` 
            - `bpcells_filtered/[sample]`: Filtered fragments with only cells passing QC
            - `bpcells(_filtered)/merged`: Fragment files containing a merge of all samples with or without QC filtering
            - `archr/[sample].arrow`, `archr/[sample].filtered.arrow`: ArchR files with and without filtering
            - `snapatac2/[sample].h5ad`, `snapatac2/filtered.h5ad`: SnapATAC2 files with and without filtering
            - `matrices/`: Peak and tile matrices in feature-major or cell-major (`_transpose`) order
                - `peaks/[sample]`, `peaks_transpose/[sample]`, `tiles/[sample]`, `tiles_transpose/[sample]`: BPCells-format compressed peak and tile matrices per-sample
                - `peaks/merged`, `peaks_transpose/merged`, `tiles/merged`, `tiles_transpose/merged`: Whole dataset (merged) peak and tile matrices
    - Task summary: 657/657 completed. Total time: 3 days, 12:31:39
- `datasets/rna-download`:
    - Download RNA datasets
    - Produces:
        - `{DATA_ROOT}/rna/downloads/[dataset]/`: Raw file downloads for each RNA benchmark dataset (excluding CELLxGENE)
    - Task summary: 9/9 completed. Total time: 3:05:25
- `datasets/rna-process`:
    - Process download RNA files into uniform inputs
    - Requires: `rna-download`
    - Produces:
        - `{DATA_ROOT}/rna/[dataset]/`:
            - `10x.h5`: 10x-format HDF5 matrix file
            - `bpcells`: Cell-major BPCells counts matrix
            - `bpcells_transpose`: Feature-major BPCells counts matrix (note the meaning of `_transpose` is flipped relative to peak and tile matrices)
            - `clusts.txt`: Cluster assignments for each cell, one cluster ID per line in same order as the cells in the data matrices. Used for marker feature benchmarks
            - `variable_genes.txt`: Pre-calculated variable gene sets to keep PCA benchmark calculations exactly matching
            - `svd.rds`: RDS file containing the dimensionality reduction results used for cluster assignments
    - Task summary: 9/9 completed. Total time: 8:14:54
- `datasets/dataset-stats`: 
    - Collect basic size statistics on each input dataset
    - Requires: `atac-process`, `rna-process`
    - Produces:
        - `{RESULTS_ROOT}/raw/datasets/dataset-stats/`: Raw outputs
        - `{REAULTS_ROOT}/data_tables/datasets/`
    - Task summary: 15/15 completed. Total time: 0:10:41