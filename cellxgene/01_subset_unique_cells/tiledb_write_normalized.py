import datetime
import gc
import os
import sys
import time

import numpy as np
import numba
import pyarrow as pa
import tiledb
import tiledbsoma

args = sys.argv[1:]
assert len(args) == 6
input_path = args[0]
output_path = args[1]
census_root = args[2]
compression_level = int(args[3])
output_stats = args[4]
buf_size_gb = float(args[5])

cpu_count = os.environ["OMP_NUM_THREADS"]
var_metadata_path = f"{census_root}/homo_sapiens/ms/RNA/var"
obs_metadata_path = f"{census_root}/homo_sapiens/obs"


def timestamp_now():
    return datetime.datetime.isoformat(datetime.datetime.now())


# Adapted from https://github.com/chanzuckerberg/cellxgene-census/blob/ec0e754b84c2ca213a1fe976c352318040c8048f/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/globals.py#L392-L406
# and from https://github.com/chanzuckerberg/cellxgene-census/blob/ec0e754b84c2ca213a1fe976c352318040c8048f/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/consolidate.py#L150
# MIT license Copyright 2022 Chan Zuckerberg Initiative (https://github.com/chanzuckerberg/cellxgene-census/blob/ec0e754b84c2ca213a1fe976c352318040c8048f/LICENSE)
context = tiledbsoma.SOMATileDBContext(
    tiledb_config={
        "py.init_buffer_bytes": int(buf_size_gb * 1024**3),
        "py.deduplicate": "true",
        "soma.init_buffer_bytes": int(buf_size_gb * 1024**3),
        "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3,
        "sm.consolidation.total_buffer_size": int(buf_size_gb * 1024**3),
        "sm.compute_concurrency_level": cpu_count,
        "sm.io_concurrency_level": cpu_count,
    }
)

platform_config = {
    "tiledb": {
        "create": {
            "capacity": 2**16,
            "dims": {
                "soma_dim_0": {
                    "tile": 2048,
                    "filters": [{"_type": "ZstdFilter", "level": compression_level}],
                },
                "soma_dim_1": {
                    "tile": 2048,
                    "filters": [
                        "ByteShuffleFilter",
                        {"_type": "ZstdFilter", "level": compression_level},
                    ],
                },
            },
            "attrs": {
                "soma_data": {
                    "filters": [
                        {
                            "_type": "ZstdFilter",
                            "level": 13
                            if compression_level == 9
                            else compression_level,
                        }
                    ]
                }
            },
            "cell_order": "row-major",
            "tile_order": "row-major",
            "allows_duplicates": True,
        },
    }
}

# Full-gene assays have special handling in the "normalized" X layers
FULL_GENE_ASSAY = [
    "EFO:0003755",  # FL-cDNA
    "EFO:0008441",  # full length single cell RNA sequencing
    "EFO:0008747",  # FRISCR
    "EFO:0008763",  # Hi-SCL
    "EFO:0008797",  # MATQ-seq
    "EFO:0008877",  # Quartz-seq
    "EFO:0008930",  # Smart-seq
    "EFO:0008931",  # Smart-seq2
    "EFO:0008956",  # SUPeR-seq
    "EFO:0009810",  # full length single nucleus RNA sequencing
    "EFO:0010004",  # SCRB-seq
    "EFO:0010022",  # Smart-3Seq
    "EFO:0010058",  # Fluidigm C1-based library preparation
    "EFO:0010184",  # Smart-like
    "EFO:0022396",  # TruSeq
    "EFO:0022488",  # Smart-seq3
    "EFO:0022839",  # STORM-seq
    "EFO:0030031",  # SCOPE-chip
    "EFO:0030061",  # mcSCRB-seq
    "EFO:0700016",  # Smart-seq v4
]


# Adapted from: https://github.com/chanzuckerberg/cellxgene-census/blob/ec0e754b84c2ca213a1fe976c352318040c8048f/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/experiment_builder.py#L774-L811
# MIT license Copyright 2022 Chan Zuckerberg Initiative (https://github.com/chanzuckerberg/cellxgene-census/blob/ec0e754b84c2ca213a1fe976c352318040c8048f/LICENSE)
@numba.jit(nopython=True, nogil=True)
def _normalize_data(
    d0,  # npt.NDArray[np.int64]
    d1,  # npt.NDArray[np.int64]
    data,  # npt.NDArray[np.float32]
    row_sum,  # npt.NDArray[np.float64]
    is_full_length,  # npt.NDArray["bool"]
    feature_length,  # np.tNDArray[np.int64]
    keepbits=15,
):
    """IMPORTANT: modifies `data` in place. Fuse the separate operators used in cellxgene_census_builder for improved performance
    1. For rows where is_full_length is true, divide by feature length
    2. Divide by row sums
    3. Reduce the data precision, with round half to even
    """
    assert data.dtype is np.dtype(np.float32)
    nmant = 23
    bits = 32
    maskbits = nmant - keepbits
    full_mask = (1 << bits) - 1
    mask = (full_mask >> maskbits) << maskbits
    half_quantum1 = (1 << (maskbits - 1)) - 1

    for i in range(len(d0)):
        x = np.float32(data[i])
        if is_full_length[d0[i]]:
            x /= feature_length[d1[i]]
        x /= row_sum[d0[i]]
        data[i] = x
    # Reduce precision
    b = data.view(np.int32)
    b += ((b >> maskbits) & 1) + half_quantum1
    b &= mask
    return data


# Back to non-CZI-derived code:
@numba.jit(nopython=True, nogil=True)
def _update_row_sum(d0, d1, data, row_sum, feature_length):
    """important: modifies `row_sum` in place"""
    for i in range(len(d0)):
        x = data[i]
        x /= feature_length[d1[i]]
        row_sum[d0[i]] += x


obs = tiledbsoma.DataFrame.open(obs_metadata_path, context=context)
obs = obs.read(
    column_names=["is_primary_data", "raw_sum", "assay_ontology_term_id"]
).concat()
feature_length = (
    tiledbsoma.DataFrame.open(var_metadata_path, context=context)
    .read(column_names=["feature_length"])
    .concat()["feature_length"]
    .to_numpy()
)
is_primary_data = obs["is_primary_data"].to_numpy()
is_full_gene_assay = np.isin(obs["assay_ontology_term_id"].to_numpy(), FULL_GENE_ASSAY)
row_sum = obs["raw_sum"].to_numpy().copy()

unique_cells = np.flatnonzero(is_primary_data)
new_idx = np.full(is_primary_data.shape, -1, dtype=np.int64)
new_idx[unique_cells] = np.arange(len(unique_cells))

input = tiledbsoma.SparseNDArray.open(
    input_path, context=context, platform_config=platform_config
)
output = tiledbsoma.SparseNDArray.create(
    output_path,
    type=pa.float32(),
    shape=(len(unique_cells), input.shape[1]),
    platform_config=platform_config,
    context=context,
)
print(f"Calculating full_gene row_sums: {timestamp_now()}", flush=True)
full_gene_data = (
    input.read(coords=[np.flatnonzero(is_primary_data & is_full_gene_assay)])
    .tables()
    .concat()
)
row_sum[is_primary_data & is_full_gene_assay] = 0
_update_row_sum(
    full_gene_data["soma_dim_0"].to_numpy(),
    full_gene_data["soma_dim_1"].to_numpy(),
    full_gene_data["soma_data"].to_numpy(),
    row_sum,
    feature_length,
)
del full_gene_data
gc.collect()


print(f"Copying data in batches: {timestamp_now()}", flush=True)
total_read = 0
elapsed_start = time.time()
process_start = time.process_time()
for batch in input.read(coords=[unique_cells]).tables():
    total_read += batch.shape[0]
    print(
        f"{timestamp_now()} Normalizing batch with size {batch.shape[0]} (total read={total_read})",
        flush=True,
    )
    new_soma_data = _normalize_data(
        batch["soma_dim_0"].to_numpy(),
        batch["soma_dim_1"].to_numpy(),
        batch["soma_data"].to_numpy().copy(),
        row_sum,
        is_full_gene_assay,
        feature_length,
    )
    new_soma_dim_0 = new_idx[batch["soma_dim_0"]]
    modified_batch = pa.Table.from_arrays(
        [new_soma_dim_0, batch["soma_dim_1"], new_soma_data],
        names=["soma_dim_0", "soma_dim_1", "soma_data"],
    )
    print(f"{timestamp_now()} Writing batch", flush=True)
    _ = output.write(modified_batch)
    print(f"{timestamp_now()} Reading next batch", flush=True)

print(f"Starting consolidate and vacuum: {timestamp_now()}", flush=True)
tiledb.consolidate(output_path)
tiledb.vacuum(output_path)
print(f"Done: {timestamp_now()}", flush=True)

cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

total_bytes = 0
for dir, _, files in os.walk(output_path):
    for f in files:
        total_bytes += os.path.getsize(os.path.join(dir, f))

# Output data
open(output_stats, "w").writelines(
    ["time_cpu\ttime_elapsed\tbytes\n", f"{cpu_time}\t{elapsed_time}\t{total_bytes}"]
)
