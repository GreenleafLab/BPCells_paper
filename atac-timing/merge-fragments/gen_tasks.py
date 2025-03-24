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


# Leave some buffer on top of the 3 hour cutoff because the unix_sort prep might take some time
THREADS=6
SBATCH_OPTS=" ".join([
    "-p biochem",
    "--constraint=CPU_SKU:7502",
    f"-c {THREADS}",
    "--time=5:00:00",
    "--mem=24G"
])

print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"

datasets = ["1m_brain", "35k_hematopoiesis", "45k_pancreas", "500k_heart"]
replicates = ["rep1","rep2","rep3","rep4","rep5"]

# List of scripts to call for each step of benchmarking for each tool

run_tools = {
    "bpcells": "Rscript merge-fragments-bpcells.R",
    "unix_sort": "python merge-fragments-unix-sort.py",
}
prep_tools = {
    "bpcells": "Rscript prep-inputs-bpcells.R",
    "unix_sort": "Rscript prep-inputs-unix-sort.R",  
}

for rep, dataset, tool in itertools.product(replicates, datasets, run_tools.keys()):
        dataset_dir = f"{DATA_ROOT}/atac/{dataset}"
        output_dir = f"{RESULTS_ROOT}/raw/atac-timing/merge-fragments/{tool}_{dataset}/"
        output_prefix = f"{output_dir}/{rep}"
        tmp_dir = f"{TMP_ROOT}/{tool}_{dataset}_{rep}"

        os.makedirs(output_dir, exist_ok=True)
        print(
            f"mkdir -p {tmp_dir}"
            f" && {SINGULARITY} {prep_tools[tool]} {dataset_dir} {tmp_dir}"
            f" && {SINGULARITY} env time -vo {output_prefix}.gnutime.txt {run_tools[tool]} {tmp_dir} {output_prefix}.tsv"
            f" && rm -r {tmp_dir}",
            file=task_txt
        )
 