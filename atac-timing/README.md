# Code to run ATAC analysis benchmarks

- `atac-timing/merge-fragments`:
    - Merge genome-sorted fragments data from multiple files (i.e. one per sample) into a single sorted file
    - Requires: `datasets/atac-process`
    - Produces:
        - `{RESULTS_ROOT}/raw/atac-timing/merge-fragments`: Raw timing outputs
        - `{RESULTS_ROOT}/data_tables/atac-timing/merge-fragments.tsv`
    - Figures: Fig S3c
    - Task summary: 40/40 completed. Total time: 23:25:45
- `atac-timing/peak-tile-timing`
    - Peak and tile matrix creation benchmark
    - Requires: `datasets/atac-process`
    - Produces:
        - `{DATA_ROOT}/atac/{dataset}/peak-subsets/`: Peak sets used for benchmarking
        - `{RESULTS_ROOT}/raw/atac-timing/peak-tile-timing/`: Raw timing outputs
        - `{RESULTS_ROOT}/data_tables/atac-timing/peak-tile-timing.tsv.gz`
    - Figures: Fig 1e-g, Fig S1 e+f
    - Task summary: 520/520 completed. Total time: 6 days, 17:06:59