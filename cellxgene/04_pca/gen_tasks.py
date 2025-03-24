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

# Server timed out at 7 hours, so adding more padding
if IS_LAPTOP:
    THREADS=16
else:
    THREADS=32
SBATCH_OPTS=" ".join([
    "-p biochem",
    "--constraint=CPU_SKU:7502",
    f"-c {THREADS}",
    "--time=10:00:00",
    "--mem=250G"
])
SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)

print(SBATCH_OPTS, file=task_txt)

CENSUS_ROOT = f"{DATA_ROOT}/cellxgene-census/2024-07-01/"

if IS_LAPTOP:
    replicates = ["rep1", "rep2"]
else:
    replicates = ["rep1", "rep2", "rep3"]

os.makedirs(f"{RESULTS_ROOT}/raw/cellxgene/04_pca/", exist_ok=True)
os.makedirs(f"{CENSUS_ROOT}/pca/", exist_ok=True)
for rep, precision, compress  in itertools.product(replicates, ["randomized", "exact"], ["compress", "uncompress"], ):
    input_path = f"{CENSUS_ROOT}/bpcells/subset_genemajor"
    
    tmp_dir = f"{TMP_ROOT}/pca_{compress}_{precision}_{rep}"
    output_dir = f"{RESULTS_ROOT}/raw/cellxgene/04_pca/"
    output_prefix = f"{output_dir}/{compress}_{precision}_{rep}"
    output_pca = f"{CENSUS_ROOT}/pca/{compress}_{precision}_{rep}.rds"

    if IS_LAPTOP:
        copy_command = f"ln -s {input_path} {tmp_dir}/full_input"
    else:
        copy_command = f"cp -r {input_path} {tmp_dir}/full_input"

    print(
        f" mkdir -p {tmp_dir} "
        f" && {copy_command} "
        f" && {SINGULARITY} env time -vo {output_prefix}.gnutime.txt "
        f"    Rscript pca.R {tmp_dir}/full_input {CENSUS_ROOT} {tmp_dir} {compress} {precision} {output_prefix}.tsv {output_pca}"
        f" && rm -r {tmp_dir}",
        file=task_txt
    )
