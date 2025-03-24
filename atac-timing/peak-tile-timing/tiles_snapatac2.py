import time
import sys

import snapatac2

args = sys.argv[1:]
assert len(args) == 4
dataset_dir = args[0]
output_timing = args[1]
staging_dir = args[2]
ref_dir = args[3]


adata = snapatac2.read(f"{staging_dir}/sample.h5ad")


elapsed_start = time.time()
process_start = time.process_time()
snapatac2.pp.add_tile_matrix(adata, n_jobs=1, exclude_chroms=["chrM"])
cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\ttool\n",
    f"{cpu_time}\t{elapsed_time}\tsnapatac2\n"
])