
# Unzip and clean up the tabula sapiens to a standard h5ad file
import sys
import subprocess

import scanpy as sc
import anndata

args = sys.argv[1:]
assert len(args) == 3
in_path = args[0]
out_path = args[1]
temp_dir = args[2]

# Unzip input file
subprocess.run(["unzip", in_path, "TabulaSapiens.h5ad", "-d", temp_dir])

# Clean up explicit zeros
d = sc.read_h5ad(f"{temp_dir}/TabulaSapiens.h5ad", backed="r")
mat = d.layers["raw_counts"][:,:]
mat.eliminate_zeros()

trimmed = anndata.AnnData(X=mat, obs = d.obs, var = d.var)
trimmed.write_h5ad(out_path)