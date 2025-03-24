import os
from pathlib import Path
import urllib.request
import re
import subprocess
file_dir = Path(__file__).parent.resolve()
task_txt = open(file_dir / "tasks.txt", "w")

root_dir = Path(__file__).parent.resolve()
while not (root_dir / "config_vars.sh").is_file() and not root_dir.parent == root_dir:
    root_dir = root_dir.parent
config = {}
exec(open(root_dir / "config_vars.sh").read(), config)

# 1.5h is usually sufficient, but leave some padding time in case of a slower node
THREADS=16
SBATCH_OPTS=" ".join([
    "-p wjg,biochem,sfgf,owners",
    f"-c {THREADS}",
    "--mem=92G",
    "--time=5:00:00",
])
print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-locked-v1", OMP_NUM_THREADS=THREADS)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]

def finish_prep(sample):
      return (
            f"mkdir -p {TMP_ROOT}/finish_prep_{sample}"
            f" && {SINGULARITY} Rscript transpose_matrix.R {DATA_ROOT}/rna/{sample}/bpcells {DATA_ROOT}/rna/{sample}/bpcells_transpose "
            f" && {SINGULARITY} Rscript prep_var_genes_clusts.R {DATA_ROOT}/rna/{sample}/bpcells {TMP_ROOT}/finish_prep_{sample} "
                f"{DATA_ROOT}/rna/{sample}/variable_genes.txt {DATA_ROOT}/rna/{sample}/svd.rds {DATA_ROOT}/rna/{sample}/clusts.txt"
            f" && rm -r {TMP_ROOT}/finish_prep_{sample}"
      )

# 1m_neurons; 7 minutes, <1GB RAM
print(f"{SINGULARITY} Rscript normalize_matrices.R {DATA_ROOT}/rna/downloads/1m_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5 "
      f" {DATA_ROOT}/rna/1m_neurons"
      f" && {finish_prep('1m_neurons')}", file=task_txt)

# 2m_perturbseq; 15 min, <16GB RAM
print(f"mkdir -p {TMP_ROOT}/2m_perturbseq "
      f" && {SINGULARITY} python 2m_perturbseq.py {DATA_ROOT}/rna/downloads/2m_perturbseq/K562_gwps_raw_singlecell_01.h5ad {TMP_ROOT}/2m_perturbseq "
      f" && {SINGULARITY} Rscript 2m_perturbseq.R {TMP_ROOT}/2m_perturbseq {DATA_ROOT}/rna/2m_perturbseq "
      f" && {finish_prep('2m_perturbseq')}"
      f" && rm -r {TMP_ROOT}/2m_perturbseq", file=task_txt)

# 4m_fetal; 12 min
print(f"{SINGULARITY} Rscript 4m_fetal.R {DATA_ROOT}/rna/downloads/4m_fetal {DATA_ROOT}/rna/4m_fetal {TMP_ROOT}"
      f" && {finish_prep('4m_fetal')}", file=task_txt)

# 11m_jax; 35m
print(f"{SINGULARITY} Rscript 11m_jax.R {DATA_ROOT}/rna/downloads/11m_jax {DATA_ROOT}/rna/11m_jax"
      f" && {finish_prep('11m_jax')}", file=task_txt)

# 20k_pbmc <1 min
print(f"{SINGULARITY} Rscript normalize_matrices.R {DATA_ROOT}/rna/downloads/20k_pbmc/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5 "
      f" {DATA_ROOT}/rna/20k_pbmc"
      f" && {finish_prep('20k_pbmc')}", file=task_txt)

# 22m_pansci 40 min
print(f" mkdir -p {DATA_ROOT}/rna/22m_pansci/"
      f" && {SINGULARITY} gunzip -c {DATA_ROOT}/rna/downloads/22m_pansci/GSE247719_20240213_PanSci_all_cells_adata.h5ad.gz > {DATA_ROOT}/rna/22m_pansci/anndata.h5ad"
      f" && {SINGULARITY} Rscript normalize_matrices.R {DATA_ROOT}/rna/22m_pansci/anndata.h5ad {DATA_ROOT}/rna/22m_pansci/ "
      f" && rm {DATA_ROOT}/rna/22m_pansci/anndata.h5ad"
      f" && {finish_prep('22m_pansci')}", file=task_txt)

# 130k_thymus_atlas 1.5 min
print(f"{SINGULARITY} python 130k_thymus_atlas.py {DATA_ROOT}/rna/downloads/130k_thymus_atlas/HTA07.A01.v02.entire_data_raw_count.h5ad {TMP_ROOT}/thymus_atlas.h5ad"
      f" && {SINGULARITY} Rscript normalize_matrices.R {TMP_ROOT}/thymus_atlas.h5ad {DATA_ROOT}/rna/130k_thymus_atlas", 
      f" && rm -r {TMP_ROOT}/thymus_atlas.h5ad"
      f" && {finish_prep('130k_thymus_atlas')}", file=task_txt)

# 480k_tabula_sapiens
print(f"mkdir -p {DATA_ROOT}/rna/480k_tabula_sapiens "
      f" && {SINGULARITY} python 480k_tabula_sapiens.py "
            f" {DATA_ROOT}/rna/downloads/480k_tabula_sapiens/TabulaSapiens.h5ad.zip"
            f" {DATA_ROOT}/rna/480k_tabula_sapiens/anndata.h5ad"
            f" {TMP_ROOT}"
      f" && {SINGULARITY} Rscript normalize_matrices.R {DATA_ROOT}/rna/480k_tabula_sapiens/anndata.h5ad {DATA_ROOT}/rna/480k_tabula_sapiens/ "
      f" && {finish_prep('480k_tabula_sapiens')}"
      f" && rm {DATA_ROOT}/rna/480k_tabula_sapiens/anndata.h5ad", file=task_txt)

# 500k_drugscreen 3 min
print(f"{SINGULARITY} Rscript normalize_matrices.R {DATA_ROOT}/rna/downloads/500k_drugscreen/H1975_A549_DrugScreen_3p_HT_nextgem_count_filtered_feature_bc_matrix.h5"
      f"  {DATA_ROOT}/rna/500k_drugscreen"
      f" && {finish_prep('500k_drugscreen')}", file=task_txt)