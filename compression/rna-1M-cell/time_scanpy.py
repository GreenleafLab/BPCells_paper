import gc
import os.path
import sys
import time

import pandas as pd
import numpy as np

import anndata as ad
import scanpy as sc

# Disable blosc multithreading for zarr for fair comparison
from numcodecs import blosc
blosc.use_threads = False


path_10x = sys.argv[1]
scratch_dir = sys.argv[2]
output_path = sys.argv[3]

assert len(sys.argv) == 4

assert path_10x.endswith(".h5")
assert os.path.isdir(scratch_dir)
assert output_path.endswith(".tsv")

def dump_mat_raw(m, p):
    f = open(p, "wb")
    m.X.data.tofile(f)
    m.X.indices.tofile(f)

methods = {
    "h5ad": {
        "path": os.path.join(scratch_dir, "mat.h5ad"),
        "read": lambda p: ad.read_h5ad(p),
        "write": lambda m, p: m.write_h5ad(p)
    },
    # Disabling this test because it's not competitive with zarr
    # "h5ad_lzf": {
    #     "path": os.path.join(scratch_dir, "mat_lzf.h5ad"),
    #     "read": lambda p: ad.read_h5ad(p),
    #     "write": lambda m, p: m.write_h5ad(p, compression='lzf')
    # },
    "h5ad_zarr": {
        "path": os.path.join(scratch_dir, "mat.zarr"),
        "read": lambda p: ad.read_zarr(p),
        "write": lambda m, p: m.write_zarr(p)
    },
    "10x": {
        "path": path_10x,
        "read": lambda p: sc.read_10x_h5(p),
        "write": None
    },
    # Disabling this test and just using later cp operations to judge I/O perf
    # "uncompressed": {
    #     "path": os.path.join(scratch_dir, "uncompressed"),
    #     "read": lambda p: open(p, "rb").read(),
    #     "write": dump_mat_raw
    # },
}

# Adapted from: https://stackoverflow.com/a/1392549
def get_size(path):
    if os.path.isfile(path):
        return os.path.getsize(path)
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)
    return total_size


data = sc.read_10x_h5(path_10x)
results = []
for name, method in methods.items():
    # Time write
    time_write = np.nan
    if method["write"] != None:
        elapsed_start = time.time()
        process_start = time.process_time()
        method["write"](data, method["path"])
        write_cpu = time.process_time() - process_start
        write_elapsed = time.time() - elapsed_start

    size = get_size(method["path"])
    # Time read
    elapsed_start = time.time()
    process_start = time.process_time()
    m2 = method["read"](method["path"])
    read_cpu = time.process_time() - process_start
    read_elapsed = time.time() - elapsed_start
    m2 = None
    gc.collect()
    results.append({
        "format": name,
        "tool": "python",
        "read": read_elapsed,
        "read_seconds_cpu": read_cpu,
        "write": write_elapsed,
        "write_seconds_cpu": write_cpu, 
        "size": size,
    })
    
df = pd.DataFrame(results)
df.to_csv(output_path, index=False, sep="\t")
