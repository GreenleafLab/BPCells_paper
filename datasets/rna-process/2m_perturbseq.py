import datetime
import gc
import sys
import scipy.sparse

import scanpy as sc
import anndata

assert len(sys.argv) == 3
in_path = sys.argv[1]
out_dir = sys.argv[2]

# Convert from the dense matrix format to a sparse format

d = sc.read_h5ad(in_path, backed='r')
chunk_size = 100000
chunk_id = 0
for i in range(0, d.X.shape[0], chunk_size):
    print(f"Loading chunk at offset {i}: {datetime.datetime.isoformat(datetime.datetime.now())}")
    sparse_mat = scipy.sparse.csr_matrix(d.X[i:i+chunk_size,:])
    adata = anndata.AnnData(X=sparse_mat, obs=d.obs.iloc[i:i+chunk_size,:], var=d.var)
    adata.write_h5ad(f"{out_dir}/{chunk_id:02d}.h5ad")
    chunk_id += 1
    gc.collect()
