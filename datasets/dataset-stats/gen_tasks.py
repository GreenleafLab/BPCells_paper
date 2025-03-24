import itertools
import os
from pathlib import Path

file_dir = Path(__file__).parent.resolve()
task_txt = open(file_dir / "tasks.txt", "w")

root_dir = Path(__file__).parent.resolve()
while not (root_dir / "config_vars.sh").is_file() and not root_dir.parent == root_dir:
    root_dir = root_dir.parent
config = {}
exec(open(root_dir / "config_vars.sh").read(), config)


THREADS=8
SBATCH_OPTS=" ".join([
    "-p biochem,owners,sfgf,wjg",
    f"-c {THREADS}",
    "--time=3:00:00",
    "--mem=32000M"
])

print(SBATCH_OPTS, file=task_txt)

DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)

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

os.makedirs(f"{RESULTS_ROOT}/raw/datasets/dataset-stats", exist_ok=True)
for dataset in rna_datasets:
    print(f"{SINGULARITY} Rscript rna_dataset_stats.R {DATA_ROOT}/rna/{dataset}/bpcells {RESULTS_ROOT}/raw/datasets/dataset-stats/rna-{dataset}.tsv",
          file=task_txt)

atac_datasets = ["1m_brain", "3k_pbmc", "10k_pbmc", "35k_hematopoiesis", "45k_pancreas", "500k_heart"]

os.makedirs(f"{RESULTS_ROOT}/raw/datasets/dataset-stats", exist_ok=True)
for dataset in atac_datasets:
    print(f"{SINGULARITY} Rscript atac_dataset_stats.R {DATA_ROOT}/atac/{dataset}/ {RESULTS_ROOT}/raw/datasets/dataset-stats/atac-{dataset}.tsv",
          file=task_txt)
