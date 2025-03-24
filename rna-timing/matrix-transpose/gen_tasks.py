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

# Ask for the maximum memory possible, since scipy is likely to need all it can get
THREADS=16
SBATCH_OPTS=" ".join([
    "-p biochem,owners",
    "--constraint=CPU_SKU:7502",
    "--constraint=CLASS:SH3_CBASE",
    f"-c {THREADS}",
    "--time=1:00:00",
    "--mem=256000M"
])


print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"

datasets = [
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
replicates = ["rep1","rep2","rep3","rep4", "rep5"]
tools = ["bpcells", "dgCMatrix", "scipy"]

# List of scripts to call for each step of benchmarking for each tool
run_tools = {
    "bpcells": "Rscript transpose_bpcells.R",
    "dgCMatrix": "Rscript transpose_dgCMatrix.R",
    "scipy": "python transpose_scipy.py",
}

for rep, dataset, tool in itertools.product(replicates, datasets, run_tools.keys()):
        dataset_dir = f"{DATA_ROOT}/rna/{dataset}"
        output_dir = f"{RESULTS_ROOT}/raw/rna-timing/matrix-transpose/{tool}_{dataset}/"
        output_prefix = f"{output_dir}/{rep}"
        tmp_dir = f"{TMP_ROOT}/transpose_{tool}_{dataset}_{rep}"

        os.makedirs(output_dir, exist_ok=True)
        print(
            f"mkdir -p {tmp_dir}"
            f" && {SINGULARITY} Rscript prep_inputs.R {dataset_dir} {tmp_dir} {tool}"
            f" && {SINGULARITY} env time -vo {output_prefix}.gnutime.txt {run_tools[tool]} {tmp_dir} {output_prefix}.tsv"
            f" && rm -r {tmp_dir}",
            file=task_txt
        )
 