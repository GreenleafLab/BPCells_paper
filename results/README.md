## Plotting results

Plotting scripts are present in `plots/atac-timing.R`, `plots/bitwidth-stats.R`, `plots/cellxgene.R`, `plots/compression.R`, and `plots/rna-timing.R`. To run these scripts yourself, clone this repository, then replace this line at the top of each script with the root folder of the repository:

```r
repo_root <- "/home/bparks/Sync/BPCells_paper/final_upload_version/"
```

Each script should be able to run end-to-end to produce a folder of svg plots. Within each script there are also commented lines to reproduce numerical claims from the paper or used as annotations on top of the plots.

These scripts only read data from the files in `data_tables`.

## Results data files

Results files are mostly in TSV format or gzipped TSV format. The mapping between figure panels and data table paths is as follows:

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

The per-file column structure is as follows:

(click to expand)

<details>
<summary><code>data_tables/datasets</code></summary>

- `atac.tsv`: ATAC dataset + sample stats. One row per sample, with rows where `sample == merged` covering the whole dataset.
    - `dataset`: Dataset ID
    - `sample`: Sample ID
    - `filtered`: true/false if these stats are taken pre- or post-filtering
    - `barcodes`: Number of unique barcodes (i.e. cell count when `filtered==TRUE`)
    - `total_fragments`: Total number of fragments
    - `median_fragments`: Median number of fragments per-barcode
- `rna.tsv`: RNA dataset stats. One row per dataset
    - `dataset`: Dataset ID
    - `cells`: Number of cells
    - `genes`: Number of genes
    - `median_reads`: Median reads (generally UMIs) per cell
    - `median_genes`: Median detected genes per cell
    - `nonzero_entries`: Total number of nonzero entries in the matrix (total detected genes across all cells)
    - `fraction_nonzero`: Fraction of counts matrix entries which are non-zero
    - `mean_nonzero_val`: Mean value when excluding matrix entries that are zero

</details>
<details>

<summary><code>data_tables/rna-timing/marker-genes</code></summary>

- `performance.tsv`
    - `tool`: Tool name
    - `dataset`: Dataset ID
    - `replicate`: Experiment replicate number
    - `time_cpu`: Measured CPU time (seconds)
    - `time_elapsed`: Measured elapsed time (seconds)
    - `genes`: Number of genes in the dataset
    - `clusts`: Number of clusters used for marker calculation
    - `max_rss`: Maximum memory usage (RSS) in bytes
- `accuracy.tsv`
    - `tool1`: Tool name (first tool in comparison)
    - `tool2`: Tool name (second tool in comparison)
    - `cor`: Correlation of the log10 p-values (cut off at -290)
    - `mean_abs_diff`: Mean absolute difference of the log10 p-values (cut off at -290)
    - `dataset`: Dataset ID
</details>
<details>
<summary><code>data_tables/rna-timing/matrix-transpose.tsv</code></summary>

- `tool`: Tool name
- `dataset`: Dataset ID
- `replicate`: Experiment replicate number
- `time_cpu`: Measured CPU time (seconds)
- `time_elapsed`: Measured elapsed time (seconds)
- `max_rss`: Maximum memory usage (RSS) in bytes
</details>
<details>
<summary><code>data_tables/rna-timing/pca-benchmark</code></summary>

- `performance.tsv`
    - `tool`: Tool name
    - `dataset`: Dataset ID
    - `step`: "normalize" or "pca"
    - `threads`: Number of threads/cores enabled
    - `replicate`: Experiment replicate number
    - `time_cpu`: Measured CPU time (seconds)
    - `time_elapsed`: Measured elapsed time (seconds)
    - `max_rss`: Maximum memory usage (RSS) in bytes
    - `n_ops`: Number of reported matrix-vector multiplies during PCA, if available
- `accuracy.tsv.gz`
    - `ref`: Tool for reference PCA coordinates
    - `alt`: Tool with comparison PCA coordinates
    - `dataset`: Dataset ID
    - `replicate`: Experiment replicate number
    - `PC_ref`: Reference principal component
    - `PC_alt`: Comparison principal component
    - `pearson`: Pearson correlation of the component
    - `axis`: "gene" or "cell" depending on whether the component was correlated across cells or across genes (ie embedding vs loading)
- `staged_file_size.tsv`
    - `tool`: Tool name
    - `threads`: Number of threads/cores enabled
    - `replicate`: Experiment replicate number
    - `file_bytes`: Total bytes of the staged matrix files
    - `subset_nonzeros`: Number of non-zero entries in the gene subset matrix being used for PCA
</details>
<details>
<summary><code>data_tables/atac-timing/merge-fragments.tsv</code></summary>

- `tool`: Tool name
- `dataset`: Dataset ID
- `replicate`: Experiment replicate number
- `time_cpu`: Measured CPU time (seconds)
- `time_elapsed`: Measured elapsed time (seconds)
- `max_rss`: Maximum memory usage (RSS) in bytes
</details>
<details>
<summary><code>data_tables/atac-timing/peak-tile-timing.tsv.gz</code></summary>

- `dataset`: Dataset ID
- `sample`: Sample ID within the dataset
- `tool`: Tool name
- `region_type`: tile, peaks-10, peas-1000, or peaks-100000 (number of peaks in peak set)
- `replicate`: Experiment replicate number
- `max_rss`: Maximum memory usage (RSS) in bytes
- `time_cpu`: Measured CPU time (seconds)
- `time_elapsed`: Measured elapsed time (seconds)
</details>
<details>
<summary><code>data_tables/bitwidth-stats</code></summary>

- `fragment_bitwidth.tsv.gz`: Statisitcs on bits-per-value for fragment files
    - `dataset`: Dataset ID
    - `sample`: Sample ID within the dataset
    - `filtered`: true/false if these stats are taken pre- or post-filtering
    - `field`: data field (cell, start, or end)
    - `bitwidth`: number of bits used per value
    - `chunk_count`: number of 128-value chunks with the specified bitwidth
- `fragment_storage.tsv`: Statistics on bytes used per field for fragment files
    - `dataset`: Dataset ID
    - `sample`: Sample ID within the dataset
    - `filtered`: true/false if these stats are taken pre- or post-filtering
    - `field`: data field (cell, start, end, end_max, or metadata)
    - `bytes`: number of bytes for the data field
    - `elements`: number of values stored
- `matrix_bitwidth.tsv.gz`: Statistics on bits-per-value for RNA or ATAC matrices
    - `dataset`: Dataset ID
    - `matrix`: Matrix source (peaks, tiles, or rna)
    - `cell_major`: true/false if matrix hass cell-major ordering
    - `ordering`: ordering/sorting of non-major axis (original, mean, nonzero, shuffle)
    - `field`: data field (val, index)
    - `bitwidth`: number of bits used per value
    - `chunk_count`: number of 128-value chunks with the specified bitwidth
- `matrix_storage.tsv`: Statistics on bytes used per field for matrices
    - `dataset`: Dataset ID
    - `matrix`: Matrix source (peaks, tiles, or rna)
    - `cell_major`: true/false if matrix hass cell-major ordering
    - `ordering`: ordering/sorting of non-major axis (original, mean, nonzero, shuffle)
    - `field`: data field (val, index, metadata, idxptr)
    - `bytes`: number of bytes for the data field
    - `elements`: number of values stored
- `matrix_shape.tsv`: Row and column counts for each matrix
    - `dataset`: Dataset ID
    - `matrix`: Matrix source (peaks, tiles, or rna)
    - `rows`: Number of rows
    - `cols`: Number of columns
- `matrix_value_histogram.tsv.gz`: Histogram of matrix value counts up to 1024
    - `dataset`: Dataset ID
    - `matrix`: Matrix source (peaks, tiles, or rna)
    - `value`: Counts value (1025 stands for all values >= 1025)
    - `count`: Number of matrix entries equal to `value`
</details>
<details>
<summary><code>data_tables/compression/fragments-read-write/</code></summary>

- `import.tsv.gz`
    - `dataset`: Dataset ID
    - `sample`: Sample ID within the dataset
    - `tool`: Tool name
    - `replicate`: Experiment replicate number
    - `max_rss`: Maximum memory usage (RSS) in bytes
    - `gnutime_elapsed`: Total process elapsed time
    - `time_cpu`: Measured CPU time (seconds)
    - `time_elapsed`: Measured elapsed time (seconds)
    - `bytes`: Size of imported file
- `read.tsv.gz`
    - `dataset`: Dataset ID
    - `sample`: Sample ID within the dataset
    - `tool`: Tool name
    - `replicate`: Experiment replicate number
    - `max_rss`: Maximum memory usage (RSS) in bytes
    - `gnutime_elapsed`: Total process elapsed time
    - `time_cpu`: Measured CPU time (seconds)
    - `time_elapsed`: Measured elapsed time (seconds)
    - `total_fragments`: Number of fragments read
</details>
<details>
<summary><code>data_tables/compression/in-memory-compression.tsv.gz</code></summary>

- `sample`: Sample ID
- `field`: Data field (cell/start/end for ATAC fragments, val/index for RNA counts matrix)
- `codec`: Compression algorithm name
- `level`: Compression level setting (or NA if not applicable)
- `filter`: Blosc pre-processing filter (or NA if not applicable)
- `runner`: Tool used to run benchmark (bpcells, lzbench, pyblosc2, or numpy)
- `replicate`: Experiment replicate number
- `write_min`: Minimum measured write (compression) time
- `read_min`: Minimum measured read (decompression) time
- `write_median`: Median measured write (compression) time (or NA if not recorded)
- `read_median`: Median measured read (decompression) time (or NA if not recorded)
- `write_iterations`: Number of compression trials performed
- `read_iterations`: Number of decompression trials performed
- `bytes`: Size of compressed data in bytes
- `input_bytes`: Size of uncompressed data in bytes
- `software_version`: Version number of compression software (or NA)
</details>
<details>
<summary><code>data_tables/compression/rna-1M-cell.tsv</code></summary>

- `replicate`: Experiment replicate number
- `tool`: Tool used to run benchmark (bpcells, python, or cp for disk-speed controls)
- `format`: File format
- `read`: Elapsed read time (seconds)
- `read_seconds_cpu`: CPU read time (seconds)
- `write`: Elapsed write time (seconds)
- `write_seconds_cpu`: CPU write time (seconds)
- `size`: File size in bytes
</details>
<details>
<summary><code>data_tables/cellxgene/01_subset_unique_cells-laptop.tsv</code></summary>

- `tool`: "bpcells" or "tiledb"
- `format`: Which matrix output is being produced
- `compression_level`: For TileDB, zstd compression level
- `replicate`: Experiment replicate number
- `time_cpu`: Measured CPU time (seconds)
- `time_elapsed`: Measured elapsed time (seconds)
- `bytes`: Total bytes to store the matrix layer
- `max_rss`: Maximum memory usage (RSS) in bytes
</details>
<details>
<summary><code>data_tables/cellxgene/02_matrix_slicing.tsv</code></summary>

- `tool`: "bpcells" or "tiledb"
- `format`: Which matrix format and/or normalization is being used
- `replicate`: Experiment replicate number
- `axis`: Matrix axis that is being subset -- "gene" or "cell"
- `slice_size`: Number of indices being selected
- `slice_type`: "random" or "sequential"
- `slice_id`: Which set of test indices is being used (assigned per `axis`, `slice_size` combination)
- `time_cpu`: Measured CPU time (seconds)
- `time_elapsed`: Measured elapsed time (seconds)
- `entries_loaded`: Number of non-zero entries loaded in the query
</details>
<details>
<summary><code>data_tables/cellxgene/03_mean_variance.tsv</code></summary>

- `tool`: "bpcells" or "tiledb"
- `format`: Which matrix format and/or normalization is being used
- `normalization`: What normalization type is being performed prior to statistics calculations
- `threads`: How many threads were used
- `replicate`: Experiment replicate number
- `time_cpu`: Measured CPU time (seconds)
- `time_elapsed`: Measured elapsed time (seconds)
- `axis`: Matrix axis that is being having statistics calculated -- "gene" or "cell"
</details>
<details>
<summary><code>data_tables/cellxgene/04_pca</code></summary>

- `timing.tsv`: Timing results
    - `input`: Input BPCells matrix compression type ("compress" or "uncompress")
    - `method`: "exact" (regular PCA) or "randomized"
    - `time_cpu`: Measured CPU time (seconds)
    - `time_elapsed`: Measured elapsed time (seconds)
    - `max_rss`: Maximum memory usage (RSS) in bytes
- `pca-cor.tsv.gz`: Accuracy results
    - `PC_ref`: Reference principal component
    - `PC_alt`: Comparison principal component
    - `pearson`: Pearson correlation of the component
    - `ref`: Source for reference PCA coordinates
    - `alt`: Source for comparison PCA coordinates
    - `axis`: "gene" or "cell" depending on whether the component was correlated across cells or across genes (ie embedding vs loading)
    - `axis_len`: Length of the axis used for correlation (i.e. number of genes or cells)