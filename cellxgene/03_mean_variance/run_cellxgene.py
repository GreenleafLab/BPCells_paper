import os
import time
import sys

import tiledbsoma
import cellxgene_census.experimental.pp

args = sys.argv[1:]
assert len(args) == 3
input_dir = args[0]
output_dir = args[1]
threads = args[2]

output_timing = f"{output_dir}/timing.tsv"

# Adapted from https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/consolidate.py#L14
# and from https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/globals.py#L392
context = tiledbsoma.SOMATileDBContext(tiledb_config = {
    "py.init_buffer_bytes": 1 * 1024**3,
    "py.deduplicate": "true",
    "soma.init_buffer_bytes": 1 * 1024**3,
    "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3,
    "sm.consolidation.total_buffer_size": 1 * 1024**3,
    "sm.compute_concurrency_level": threads,
    "sm.io_concurrency_level": threads,
})

census = tiledbsoma.Experiment.open(f"{input_dir}/census", context=context)

nonzeros = 0

elapsed_start = time.time()
process_start = time.process_time()
with tiledbsoma.ExperimentAxisQuery(census, measurement_name="RNA") as query:
    gene_res = cellxgene_census.experimental.pp.mean_variance(
        query, 
        layer = "layer",
        axis = 0,
        calculate_mean=True, 
        calculate_variance=True,
    )
cpu_time_gene = time.process_time() - process_start
elapsed_time_gene = time.time() - elapsed_start

elapsed_start = time.time()
process_start = time.process_time()
with tiledbsoma.ExperimentAxisQuery(census, measurement_name="RNA") as query:
    cell_res = cellxgene_census.experimental.pp.mean_variance(
        query, 
        layer = "layer",
        axis = 1,
        calculate_mean=True, 
        calculate_variance=True,
    )
cpu_time_cell = time.process_time() - process_start
elapsed_time_cell = time.time() - elapsed_start

# Output data
open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\taxis\n",
    f"{cpu_time_cell}\t{elapsed_time_cell}\tcell\n",
    f"{cpu_time_gene}\t{elapsed_time_gene}\tgene",
])

gene_res.to_csv(f"{output_dir}/gene_stats.tsv", sep="\t", index=False)
cell_res.to_csv(f"{output_dir}/cell_stats.tsv", sep="\t", index=False)