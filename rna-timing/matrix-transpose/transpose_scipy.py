import gc

import sys
import time


import h5py
import scipy
import numpy as np

assert len(sys.argv) == 3

tmp_dir = sys.argv[1]
output_timing = sys.argv[2]

with h5py.File(f"{tmp_dir}/input.h5") as f:
    mat = scipy.sparse.csc_matrix(
        (f["matrix/data"][:],  f["matrix/indices"].astype(np.int32)[:],  f["matrix/indptr"][:]),
        shape=f["matrix/shape"][:]
    )

gc.collect()

elapsed_start = time.time()
process_start = time.process_time()
new_mat = mat.tocsr()
cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\n",
    f"{cpu_time}\t{elapsed_time}\n"
])
