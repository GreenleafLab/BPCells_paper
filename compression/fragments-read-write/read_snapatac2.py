from pathlib import Path
import sys
import time

import snapatac2

args = sys.argv[1:]
assert len(args) == 4
dataset_dir = Path(args[0])
reference_dir = Path(args[1])
staging_dir = Path(args[2])
output_timing = args[3]

elapsed_start = time.time()
process_start = time.process_time()
frags = snapatac2.read(staging_dir / "sample.h5ad").obsm["fragment_paired"]
cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

fragment_count = frags.nnz

open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\ttotal_fragments\n",
    f"{cpu_time}\t{elapsed_time}\t{fragment_count}\n"
])