

# Data download
Running the CELLxGENE census dataset requires downloading a large amount of data, which will potentially incur large costs to CZI who pays to host the data on S3. These benchmarks require downloading a bit under 300GB of data, which at Amazon's listed data transfer prices would cost CZI $15-$27. Please be respectful of CZI's resources by avoiding repeated data downloads when possible.

Because of the potential costs to CZI of data downloads, we do not provide an automated script to perform the download. Please run the following commands manually to perform the download (starting from your directory containing `config_vars.sh`):

```sh
source config_vars.sh
# Check that $DATA_ROOT exists
if [ -d "${DATA_ROOT}" ]; then
    mkdir -p "$DATA_ROOT/cellxgene-census/2024-07-01"
    # The big data download:
    aws s3 sync \
        --exclude '*/RNA/X/normalized/*' \
        s3://cellxgene-data-public/cell-census/2024-07-01/soma/census_data/homo_sapiens \
        "$DATA_ROOT/cellxgene-census/2024-07-01/homo_sapiens"

    # Smaller metadata download (runs fast)
    aws s3 sync s3://cellxgene-data-public/cell-census/2024-07-01/soma/census_info "$DATA_ROOT/cellxgene-census/2024-07-01/census_info"
fi
```

# Job descriptions

- `cellxgene/README.md` Contains instructions for downloading the CELLxGENE census dataset
    - Download Produces:
        - `{DATA_ROOT}/cellxgene-census/2024-07-01/homo_sapiens`: Copy of the raw counts data and obs/var metadata (~300GB, mainly for the counts matrix)
        - `{DATA_ROOT}/cellxgene-census/2024-07-01/census_info`: Per-dataset metadata (relatively small)
- `cellxgene/00_convert_bpcells_chunks`:
    - Convert the TileDB input files to formats that are directly readable with BPCells
    - Requires: `cellxgene/README.md` data download complete
    - Produces:
        - `{DATA_ROOT}/cellxgene-census/2024-07-01/bpcells/cellmajor_chunks`: Copy of the full CELLxGENE census in BPCells matrix format in chunks of 100K cells
        - `{DATA_ROOT}/cellxgene-census/2024-07-01/`
            - `obs.rds`: R conversion of census dataframe `homo_sapiens/ms/RNA/obs`
            - `obs_idx.rds`: List of 1-based indices for cells where `is_primary_data` is True (i.e. R indices of the unique cells)
            - `obs_full_gene_mask.rds`: Boolean vector length # unique cells marking which cells were collected with assays that require
                gene length correction during normalization (same ordering as `obs_idx.rds`)
            - `var.rds`: R conversion of census dataframe `homo_sapiens/ms/RNA/var`
            - `feature_presence.rds`: R conversion of census matrix `homo_sapiens/ms/RNA/feature_dataset_presence_matrix`
    - Task summary:
        - Server: 76/76 completed. Total time: 23:17:46
        - Laptop: 76/76 completed. Total time: 10:50:55
- `cellxgene/01_subset_unique_cells`:
    - Run benchmarking of writing the census matrix in BPCells and TileDB formats
    - Requires: `00_convert_bpcells_chunks`
    - Produces:
        - `{DATA_ROOT}/cellxgene-census/2024-07-01/`
            - `bpcells`
                - `subset_cellmajor`: Unique cell raw counts matrix in BPCells format with cell-major ordering and 100K cell chunks
                - `subset_genemajor`: Unique cell raw counts matrix in BPCells format with gene-major ordering and 100K cell chunks
                - `subset_norm`: An RDS file containing the normalization factors for each cell
            - `tiledb`
                - `raw_zstd_1`, `raw_zstd_9`: Unique cell raw counts matrix in TileDB format, with zstd level 9 compression (default), or level 1 as a performance test.
                - `normalized_zstd_1`, `normalized_zstd_9`: Unique cell normalized counts matrix in TileDB format, with zstd level 9 (default) or level 1 compression
        - `{RESULTS_ROOT}/raw/cellxgene/01_subset_unique_cells`: Raw timing outputs
        - `{RESULTS_ROOT}/data_tables/cellxgene/01_subset_unique_cells-laptop.tsv`
    - Figures: Fig 3a+b, Fig S4a
    - Task summary:
        - Server: 21/21 completed. Total time: 2 days, 9:09:30
        - Laptop: 10/10 completed. Total time: 1 day, 15:50:04
- `cellxgene/02_matrix_slicing`:
    - Test speed of random and sequential subset queries by gene or by cell
    - Requires: `01_subset_unique_cells`
    - Produces:
        - `{DATA_ROOT}/cellxgene-census/2024-07-01/slice_coords/`: List of coordinates used for subsetting
        - `{RESULTS_ROOT}/raw/cellxgene/02_matrix_slicing`: Raw timing outputs
        - `{RESULTS_ROOT}/data_tables/cellxgene/02_matrix_slicing.tsv`
    - Figures: Fig 3c+d, Fig S4 b,e-g
    - Task summary:
        - Server: 18/18 completed. Total time: 9:10:42
        - Laptop: 12/12 completed. Total time: 8:22:36
- `cellxgene/03_mean_variance`:
    - Measure time to calculate per-gene mean and variance of the unique cells matrix
    - Requires: `01_subset_unique_cells`
    - Produces:
        - `{RESULTS_ROOT}/raw/cellxgene/03_mean_variance`: Raw timing outputs and calculated cell/gene statistics
        - `{RESULTS_ROOT}/data_tables/cellxgene/03_mean_variance.tsv`
    - Figures: Fig 3e, Fig S4 c+d
    - Task summary:
        - Server: 24/24 completed. Total time: 6:53:12
        - Laptop: 16/16 completed. Total time: 3:55:53
- `cellxgene/04_pca`:
    - Measure time to run PCA on the unique cells matrix
    - Requires: `01_subset_unique_cells`
    - Produces:
        - `{DATA_ROOT}/cellxgene-census/2024-07-01/pca`: PCA loadings and embeddings in rds format
        - `{RESULTS_ROOT}/raw/cellxgene/04_pca`: Raw timing outputs
        - `{RESULTS_ROOT}/data_tables/cellxgene/04_pca`
    - Figures: Fig 3f-h 
    - Task summary:
        - Server: 12/12 completed. Total time: 1 day, 9:05:20
        - Laptop: 8/8 completed. Total time: 1 day, 21:11:05