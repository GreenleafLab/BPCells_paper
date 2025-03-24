import datetime
import os
import sys
import time

import numpy as np
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
obs_metadata_path = f"{census_root}/homo_sapiens/obs"

def timestamp_now():
    return datetime.datetime.isoformat(datetime.datetime.now())

# Adapted from https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/consolidate.py#L14
# and from https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/globals.py#L392
context = tiledbsoma.SOMATileDBContext(tiledb_config = {
    "py.init_buffer_bytes": int(buf_size_gb * 1024**3),
    "py.deduplicate": "true",
    "soma.init_buffer_bytes": int(buf_size_gb * 1024**3),
    "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3,
    "sm.consolidation.total_buffer_size": int(buf_size_gb * 1024**3),
    "sm.compute_concurrency_level": cpu_count,
    "sm.io_concurrency_level": cpu_count,
})


platform_config =  {
    "tiledb": {
        "create": {
            "capacity": 2**16,
            "dims": {
                "soma_dim_0": {"tile": 2048, "filters": [{"_type": "ZstdFilter", "level": compression_level}]},
                "soma_dim_1": {"tile": 2048, "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": compression_level}]},
            },
            "attrs": {"soma_data": {"filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": compression_level}]}},
            "cell_order": "row-major",
            "tile_order": "row-major",
            "allows_duplicates": True,
        },
    }
}

obs = tiledbsoma.DataFrame.open(obs_metadata_path, context=context)
obs = obs.read(column_names=["is_primary_data", "raw_sum"]).concat()
is_primary_data = obs["is_primary_data"].to_numpy()
unique_cells = np.flatnonzero(is_primary_data)
new_idx = np.full(is_primary_data.shape, -1, dtype=np.int64)
new_idx[unique_cells] = np.arange(len(unique_cells))


input = tiledbsoma.SparseNDArray.open(input_path, context=context, platform_config=platform_config)

output = tiledbsoma.SparseNDArray.create(
    output_path,
    type = pa.float32(),
    shape = (len(unique_cells), input.shape[1]),
    platform_config=platform_config,
    context=context,
)

print(f"Copying data in batches: {timestamp_now()}", flush=True)
total_read = 0
elapsed_start = time.time()
process_start = time.process_time()
for batch in input.read(coords=[unique_cells]).tables():
    total_read += batch.shape[0]
    print(f"{timestamp_now()} Modifying batch with size {batch.shape[0]} (total read={total_read})", flush=True)
    new_soma_dim_0 = new_idx[batch["soma_dim_0"]]
    modified_batch = pa.Table.from_arrays(
        [new_soma_dim_0, batch["soma_dim_1"], batch["soma_data"]],
        names=["soma_dim_0", "soma_dim_1", "soma_data"]
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
open(output_stats, "w").writelines([
    "time_cpu\ttime_elapsed\tbytes\n",
    f"{cpu_time}\t{elapsed_time}\t{total_bytes}"
])
