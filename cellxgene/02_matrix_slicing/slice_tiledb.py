import time
import sys

import tiledbsoma

args = sys.argv[1:]
assert len(args) == 5
input_dir = args[0]
coords = args[1]
axis = args[2]
output_timing = args[3]
threads = int(args[4])

assert axis in ["cell", "gene"]

# Adapted from https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/consolidate.py#L14
# and from https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/globals.py#L392
context = tiledbsoma.SOMATileDBContext(tiledb_config = {
    "py.init_buffer_bytes": 4 * 1024**3,
    "py.deduplicate": "true",
    "soma.init_buffer_bytes": 4 * 1024**3,
    "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3,
    "sm.consolidation.total_buffer_size": 4 * 1024**3,
    "sm.compute_concurrency_level": threads,
    "sm.io_concurrency_level": threads,
})

mat = tiledbsoma.SparseNDArray.open(f"{input_dir}/matrix", context=context)

coords = [int(x.strip())-1 for x in open(coords, 'r')]

if axis == "cell" and len(coords) < mat.shape[0]:
    print("Subsetting on cell axis")
    coords = [coords]
elif axis == "gene" and len(coords) < mat.shape[1]:
    print("Subsetting on gene axis")
    coords = [slice(None), coords]
else:
    print("No subsetting")
    coords = []

nonzeros = 0
elapsed_start = time.time()
process_start = time.process_time()
for batch in mat.read(coords=coords).tables():
    nonzeros += len(batch)

cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

# Output data
open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\tentries_loaded\n",
    f"{cpu_time}\t{elapsed_time}\t{nonzeros}"
])