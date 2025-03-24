import os
import os.path
import shutil
import sys
import subprocess
import tempfile
import time
import uuid


def main():
    if len(sys.argv) != 6:
        print("Wrong number of args")
        sys.exit(1)

    sample_dir = sys.argv[1]
    data_dir = sys.argv[2]
    results_dir = sys.argv[3]
    tmp_dir = sys.argv[4]
    step = sys.argv[5]

    # Needed since singularity doesn't have the default cache dir writeable
    os.environ["NUMBA_CACHE_DIR"] = tmp_dir
    import scanpy as sc



    input_matrix = f"{sample_dir}/bpcells"
    input_var_genes = f"{sample_dir}/variable_genes.txt"
    assert step in ["normalize", "pca"]


    normalize_path = os.path.join(tmp_dir, "normalize-mat.h5ad")
    pca_path = os.path.join(data_dir, "pca-mat.h5ad")

    normalize_timing = os.path.join(results_dir, "normalize-timing.tsv")
    pca_timing = os.path.join(results_dir, "pca-timing.tsv")

    # Overview of the analysis if it weren't split up into steps with additional timing
    # adata = sc.read_10x_h5(input_matrix)
    # adata.var_names = adata.var.gene_ids
    # sc.pp.filter_genes(adata, min_counts=1)
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # sc.pp.scale(adata, zero_center=False)
    # sc.pp.pca(adata, zero_center=True, use_highly_variable=True)

    # Code from Scanpy example: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/dask.html
    import h5py
    import anndata as ad
    import dask
    import dask.array as da
    import dask.distributed as dd
    from scipy import sparse
    import numpy as np

    def read_sparse_as_dask(file_pth: str, elem_name: str, stride: int):
        with h5py.File(file_pth, "r") as f:
            elem = f[elem_name]
            shape = elem.attrs["shape"]
            if (encoding_type := elem.attrs["encoding-type"]) != "csr_matrix":
                raise ValueError(
                    f"This method was only written for csr_matrix encoding, but a {encoding_type} encoding was found."
                )
            dtype = elem["data"].dtype
        def make_dask_chunk(block_id=None):
            # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
            # https://github.com/scverse/anndata/issues/1105
            with h5py.File(file_pth, "r") as f:
                mtx = ad.experimental.sparse_dataset(f[elem_name])
                (row, _) = block_id
                chunk = mtx[
                    slice(
                        row * stride,
                        min((row * stride) + stride, shape[0]),
                    )
                ]
            return chunk
        chunks_0 = (stride,) * (shape[0] // stride)
        chunks_0 += (shape[0] % stride,)
        chunks_1 = (shape[1],)
        da_mtx = da.map_blocks(
            make_dask_chunk,
            dtype=dtype,
            chunks=(chunks_0, chunks_1),
            meta=sparse.csr_matrix((0, 0), dtype=np.float32),
        )
        return da_mtx

    SPARSE_CHUNK_SIZE = 100_000
    DENSE_CHUNK_SIZE = 10_000

    # I got some issues with workers being shut down from too much memory usage,
    # using the default "auto" settings.
    # This allows up to 4x memory oversubscription on the theory that high memory usage in a worker
    # might be transient.
    client = dd.Client(n_workers=int(os.environ["OMP_NUM_THREADS"]), 
                       threads_per_worker=1,
                       memory_limit=min(1.0, 4.0/int(os.environ["OMP_NUM_THREADS"])))
    adata_path = f"{tmp_dir}/input.h5ad"

    if step == "normalize":
        # Use BPCells to convert from 10x matrix to anndata format, but don't time this part
        # EDIT: this has been moved out of the normalize script to allow better checking of total cpu by /bin/time
        # subprocess.run(["Rscript", "-e", f"library(BPCells); open_matrix_dir(\"{input_matrix}\") |> convert_matrix_type(\"float\") |> write_matrix_anndata_hdf5(\"{adata_path}\")"], 
                    # check=True)
        with h5py.File(adata_path, "r") as f:
            adata = ad.AnnData(
                obs=ad.experimental.read_elem(f["obs"]),
                var=ad.experimental.read_elem(f["var"]),
            )
        adata.X = read_sparse_as_dask(adata_path, "X", SPARSE_CHUNK_SIZE)

            
        elapsed_start = time.time()
        process_start = time.process_time()
        sc.pp.normalize_total(adata, target_sum=1e4)
        cpu_time = time.process_time() - process_start
        elapsed_time = time.time() - elapsed_start
        # Manually set the variable genes
        variable_genes = set(x.strip() for x in open(input_var_genes))
        variable_genes = [x for x in adata.var.index if x in variable_genes]
        elapsed_start = time.time()
        process_start = time.process_time()
        # Do count time to calculate per-gene mean+variance
        dask.compute(sc.pp._utils._get_mean_var(adata.X))
        adata = adata[:,variable_genes]
        sc.pp.log1p(adata)
        sc.pp.scale(adata, zero_center=False)
        cpu_time += time.process_time() - process_start
        elapsed_time += time.time() - elapsed_start
        
        open(normalize_timing, "w").writelines([
            "time_cpu\ttime_elapsed\tstep\tn_ops\n",
            f"{cpu_time}\t{elapsed_time}\tnormalize\tNA\n"
        ])
    elif step == "pca":
        # Re-create the in-memory dask object outside of the timing block
        with h5py.File(adata_path, "r") as f:
            adata = ad.AnnData(
                obs=ad.experimental.read_elem(f["obs"]),
                var=ad.experimental.read_elem(f["var"]),
            )
        adata.X = read_sparse_as_dask(adata_path, "X", SPARSE_CHUNK_SIZE)
        sc.pp.normalize_total(adata, target_sum=1e4)
        variable_genes = set(x.strip() for x in open(input_var_genes))
        variable_genes = [x for x in adata.var.index if x in variable_genes]
        adata = adata[:,variable_genes]
        sc.pp.log1p(adata)
        sc.pp.scale(adata, zero_center=False)
        adata.layers["dense"] = adata.X.rechunk((DENSE_CHUNK_SIZE, -1)).map_blocks(
            lambda x: x.toarray(), dtype=adata.X.dtype, meta=np.array([])
        )
        
        elapsed_start = time.time()
        process_start = time.process_time()
        sc.pp.pca(adata, zero_center=True, layer="dense")
        adata.obsm["X_pca"] = adata.obsm["X_pca"].compute()
        cpu_time = time.process_time() - process_start
        elapsed_time = time.time() - elapsed_start

        adata_out = ad.AnnData(obs=adata.obs, var=adata.var, obsm=adata.obsm, varm=adata.varm)
        adata_out.write_h5ad(pca_path)
        open(pca_timing, "w").writelines([
            "time_cpu\ttime_elapsed\tstep\tn_ops\n",
            f"{cpu_time}\t{elapsed_time}\tpca\tNA\n"
        ])

if __name__ == "__main__":
    main()