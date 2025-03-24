import itertools
import math
import os
from pathlib import Path
import subprocess

file_dir = Path(__file__).parent.resolve()
task_txt = open(file_dir / "tasks.txt", "w")

root_dir = Path(__file__).parent.resolve()
while not (root_dir / "config_vars.sh").is_file() and not root_dir.parent == root_dir:
    root_dir = root_dir.parent
config = {}
exec(open(root_dir / "config_vars.sh").read(), config)


DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
IS_LAPTOP=bool(config["IS_LAPTOP"])

if IS_LAPTOP:
    THREADS=16
else:
    THREADS=32

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)

SBATCH_OPTS=" ".join([
    "-p biochem",
    "--constraint=CPU_SKU:7502",
    f"-c {THREADS}",
    "--time=2:00:00",
    "--mem=128G"
])

print(SBATCH_OPTS, file=task_txt)

CENSUS_ROOT = f"{DATA_ROOT}/cellxgene-census/2024-07-01/"

if IS_LAPTOP:
    replicates = ["rep1", "rep2"]
else:
    replicates = ["rep1", "rep2", "rep3"]

# Input format, slice command
slice_targets = [
    ("bpcells_cellmajor", "Rscript prep_inputs.R", "Rscript run_bpcells.R"),
    ("bpcells_genemajor", "Rscript prep_inputs.R", "Rscript run_bpcells.R"),
    ("bpcells_cellmajor_norm", "Rscript prep_inputs.R", "Rscript run_bpcells.R"),
    ("bpcells_genemajor_norm", "Rscript prep_inputs.R", "Rscript run_bpcells.R"),
    ("bpcells_cellmajor_lognorm", "Rscript prep_inputs.R", "Rscript run_bpcells.R"),
    ("bpcells_genemajor_lognorm", "Rscript prep_inputs.R", "Rscript run_bpcells.R"),
    ("tiledb", "python prep_inputs.py", "python run_cellxgene.py"),
    ("tiledb_norm", "python prep_inputs.py", "python run_cellxgene.py"),
]

thread_variants = [THREADS]

for threads, (tool, prep_command, run_command), rep in itertools.product(thread_variants, slice_targets, replicates):
    # Prep inputs once for each tool/replicate combination, then run all subsets on the same node
    staging_dir = f"{TMP_ROOT}/mean_variance_{tool}_{threads}_{rep}"
    result_dir = f"{RESULTS_ROOT}/raw/cellxgene/03_mean_variance/{tool}_{threads}_{rep}"
    os.makedirs(result_dir, exist_ok=True)

    print(
        f"mkdir -p {staging_dir}"
        f" && {SINGULARITY} {prep_command} {tool} {CENSUS_ROOT} {staging_dir} {IS_LAPTOP}"
        f" && {SINGULARITY} env time -vo {result_dir}/timing.gnutime.txt {run_command} {staging_dir} {result_dir} {threads}"
        f" && rm -r {staging_dir}",
        file=task_txt
    )




