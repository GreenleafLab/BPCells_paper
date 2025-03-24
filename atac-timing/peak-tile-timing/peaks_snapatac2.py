import time
import sys

import snapatac2

args = sys.argv[1:]
assert len(args) == 4
dataset_dir = args[0]
output_timing = args[1]
peaks_file = args[2]
staging_dir = args[3]

adata = snapatac2.read(f"{staging_dir}/sample.h5ad")

elapsed_start = time.time()
process_start = time.process_time()
snapatac2.pp.make_peak_matrix(adata, peak_file=peaks_file, file=f"{staging_dir}/matrix.h5ad")
cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

peak_count = sum(1 for _ in open(peaks_file))

open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\ttool\tpeaks\n",
    f"{cpu_time}\t{elapsed_time}\tsnapatac2\t{peak_count}\n"
])