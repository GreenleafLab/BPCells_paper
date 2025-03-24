import os
import os.path
import shutil
import sys
import tempfile
import time
import uuid


if len(sys.argv) != 6:
    print("Wrong number of args")
    sys.exit(1)

sample_dir = sys.argv[1]
data_dir = sys.argv[2]
results_dir = sys.argv[3]
tmp_dir = sys.argv[4]
step = sys.argv[5]

input_matrix = f"{sample_dir}/10x.h5"
input_var_genes = f"{sample_dir}/variable_genes.txt"

# Needed since singularity doesn't have the default cache dir writeable
os.environ["NUMBA_CACHE_DIR"] = tmp_dir
import scanpy as sc

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


if step == "normalize":
    tmp_path = os.path.join(tmp_dir, "input.h5")
    shutil.copyfile(input_matrix, tmp_path)
    adata = sc.read_10x_h5(tmp_path)
    
    elapsed_start = time.time()
    process_start = time.process_time()
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    cpu_time = time.process_time() - process_start
    elapsed_time = time.time() - elapsed_start
    # Manually set the variable genes
    variable_genes = set(x.strip() for x in open(input_var_genes))
    variable_genes = [x for x in adata.var.index if x in variable_genes]
    elapsed_start = time.time()
    process_start = time.process_time()
    # Do count time to calculate per-gene mean+variance
    sc.pp._utils._get_mean_var(adata.X)
    adata = adata[:,variable_genes]
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False)
    cpu_time += time.process_time() - process_start
    elapsed_time += time.time() - elapsed_start

    adata.write(normalize_path)
    
    open(normalize_timing, "w").writelines([
        "time_cpu\ttime_elapsed\tstep\tn_ops\n",
        f"{cpu_time}\t{elapsed_time}\tnormalize\tNA\n"
    ])
elif step == "pca":
    adata = sc.read(normalize_path)
    elapsed_start = time.time()
    process_start = time.process_time()
    sc.pp.pca(adata, zero_center=True)
    cpu_time = time.process_time() - process_start
    elapsed_time = time.time() - elapsed_start

    import anndata as ad
    adata_out = ad.AnnData(obs=adata.obs, var=adata.var, obsm=adata.obsm, varm=adata.varm)
    adata_out.write_h5ad(pca_path)

    open(pca_timing, "w").writelines([
        "time_cpu\ttime_elapsed\tstep\tn_ops\n",
        f"{cpu_time}\t{elapsed_time}\tpca\tNA\n"
    ])
