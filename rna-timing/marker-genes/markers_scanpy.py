import os
import shutil
import sys
import tempfile
import time
import uuid

import scanpy as sc
import numpy as np
import pandas as pd



if len(sys.argv) != 5:
    print("Wrong number of args")
    sys.exit(1)

method = sys.argv[1]
tmp_dir = sys.argv[2]
output_timing = sys.argv[3]
output_results = sys.argv[4]

tie_correct = method == "tiecorrect"

adata = sc.read_10x_h5(os.path.join(tmp_dir, "input.h5"))

adata.obs["clust"] = [line.strip() for line in open(os.path.join(tmp_dir, "clusts.txt"))]

sc.pp.normalize_total(adata)

elapsed_start = time.time()
process_start = time.process_time()
sc.tl.rank_genes_groups(adata, 'clust', method="wilcoxon", tie_correct=tie_correct)
cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

res = adata.uns["rank_genes_groups"]
df = pd.DataFrame({k:np.array(list(zip(*res[k]))).reshape(-1) for k in ["names", "scores", "pvals", "logfoldchanges"]})
df["clust"] = np.repeat(np.array(adata.obs["clust"].values.categories).reshape((-1, 1)), adata.shape[1], axis=1).reshape(-1)

df.to_csv(output_results, index=False)

clust_count = adata.obs["clust"].values.categories.size

open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\tgenes\tclusts\n",
    f"{cpu_time}\t{elapsed_time}\t{adata.shape[1]}\t{clust_count}\n"
])
