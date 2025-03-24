import itertools
import os
from pathlib import Path
import random
import subprocess

file_dir = Path(__file__).parent.resolve()
task_txt = open(file_dir / "tasks.txt", "w")

root_dir = Path(__file__).parent.resolve()
while not (root_dir / "config_vars.sh").is_file() and not root_dir.parent == root_dir:
    root_dir = root_dir.parent
config = {}
exec(open(root_dir / "config_vars.sh").read(), config)

# Note: Leave plenty of time and memory for matrix stats on the biggest matrices, but
# I expect most stuff will run a lot faster
THREADS=1
SBATCH_OPTS=" ".join([
    "-p biochem,owners,wjg",
    f"-c {THREADS}",
    "--time=1:30:00",
    "--mem=24G"
])

print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"

atac_datasets = [
    "1m_brain",
    "3k_pbmc",
    "10k_pbmc",
    "35k_hematopoiesis",
    "45k_pancreas",
    "500k_heart",
]
rna_datasets = [
    "20k_pbmc", 
    "130k_thymus_atlas", 
    "480k_tabula_sapiens",
    "500k_drugscreen", 
    "1m_neurons", 
    "2m_perturbseq", 
    "4m_fetal", 
    "11m_jax", 
    "22m_pansci", 
]

# Get list of (matrix_id, matrix_path)
atac_matrices = [
    (f"{dataset}__{mat_type}", f"{DATA_ROOT}/atac/{dataset}/matrices/{mat_type}/merged")
    for dataset in atac_datasets 
    for mat_type in ["peaks", "tiles", "peaks_transpose", "tiles_transpose"]
]

rna_matrices = [
    (f"{dataset}__rna{transpose}", f"{DATA_ROOT}/rna/{dataset}/bpcells{transpose}")
    for dataset in rna_datasets
    for transpose in ["", "_transpose"]
]


# For each matrix, calculate bitwidth, size info, value histogram
for (matrix_id, matrix_path) in rna_matrices + atac_matrices:
    tmp_dir = f"{TMP_ROOT}/bitwidth-stats-matrix-{matrix_id}"
    output_dir = f"{RESULTS_ROOT}/raw/compression/bitwidth-stats/matrices/{matrix_id}"
    os.makedirs(output_dir, exist_ok=True)
    print(
        f"mkdir -p {tmp_dir}"
        f" && {SINGULARITY} Rscript matrix_bitwidth_stats.R {matrix_path} {tmp_dir} {output_dir}"
        f" && rm -r {tmp_dir}",
        file=task_txt
    )

for dataset in atac_datasets:
    tmp_dir = f"{TMP_ROOT}/bitwidth-stats-fragments-{dataset}"
    input_dir = f"{DATA_ROOT}/atac/{dataset}/bpcells"
    output_dir = f"{RESULTS_ROOT}/raw/compression/bitwidth-stats/fragments/{dataset}"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir + "_filtered", exist_ok=True)
    print(
        f" {SINGULARITY} Rscript fragment_bitwidth_stats.R {input_dir} {output_dir}"
        f" && {SINGULARITY} Rscript fragment_bitwidth_stats.R {input_dir}_filtered {output_dir}_filtered",
        file=task_txt
    )

# Quick combined job to capture all matrix dimensions for density calculations
matrix_dimension_jobs = []
for (matrix_id, matrix_path) in rna_matrices + atac_matrices:
    output_dir = f"{RESULTS_ROOT}/raw/compression/bitwidth-stats/matrices/{matrix_id}"
    matrix_dimension_jobs.append(
        f"{SINGULARITY} Rscript matrix_dimension.R {matrix_path} {output_dir}",
    )
print(" && ".join(matrix_dimension_jobs), file=task_txt)