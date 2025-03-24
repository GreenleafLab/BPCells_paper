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
THREADS=3
SBATCH_OPTS=" ".join([
    "-p biochem,owners",
    "--constraint=CPU_SKU:7502",
    "--constraint=CLASS:SH3_CBASE",
    f"-c {THREADS}",
    "--time=3:00:00",
    "--mem=120G"
])


print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"

replicates = ["rep1","rep2","rep3","rep4", "rep5"]

input_paths = {
    "10k_pbmc": f"{DATA_ROOT}/atac/10k_pbmc/bpcells_filtered/10k_pbmc", # ATAC
    "130k_thymus": f"{DATA_ROOT}/rna/130k_thymus_atlas/bpcells",
    "130k_thymus_transpose": f"{DATA_ROOT}/rna/130k_thymus_atlas/bpcells_transpose",
}

input_type = {
    "10k_pbmc": "fragments",
    "130k_thymus": "matrix",
    "130k_thymus_transpose": "matrix",
}

fields = {
    "10k_pbmc": ["start", "end", "cell"],
    "130k_thymus": ["index", "val"],
    "130k_thymus_transpose": ["index", "val"],
}

tmp_dir_id = 0
# Jobs for lzbench
LZBENCH_CMD="/lzbench-master/lzbench -ezstd/zstd_fast/lz4/lz4fast,1,2,3,4,5,6,7,8,9/lz4hc,1,2,3,4,5,6,7,8,9/zlib -t2,2 -o4 -z"
for rep in replicates:
    output_dir = f"{RESULTS_ROOT}/raw/compression/in-memory-compression/lzbench_{rep}/"
    os.makedirs(output_dir, exist_ok=True)

    for input in input_paths.keys():
        for field in fields[input]:
            tmp_dir = f"{TMP_ROOT}/in-memory-compression_{tmp_dir_id}"
            tmp_dir_id += 1
            print(
                f"mkdir -p {tmp_dir} "
                f" && {SINGULARITY} Rscript decompress-bpcells.R {input_paths[input]} {input_type[input]} {tmp_dir}/input"
                f" && {SINGULARITY} {LZBENCH_CMD} {tmp_dir}/input/{field} > {output_dir}/{field}-{input} "
                f" && rm -r {tmp_dir}",
                file=task_txt
            )

# Jobs for blosc
blosc_filters = ["none", "byte", "bit", "delta", "delta_byte"]
for rep, filter in itertools.product(replicates, blosc_filters):
    output_dir = f"{RESULTS_ROOT}/raw/compression/in-memory-compression/blosc_{rep}/"
    os.makedirs(output_dir, exist_ok=True)

    for input in input_paths.keys():
        for field in fields[input]:
            tmp_dir = f"{TMP_ROOT}/in-memory-compression_{tmp_dir_id}"
            tmp_dir_id += 1
            print(
                f"mkdir -p {tmp_dir} "
                f" && {SINGULARITY} Rscript decompress-bpcells.R {input_paths[input]} {input_type[input]} {tmp_dir}/input"
                f" && {SINGULARITY} python compress_blosc.py "
                    f" {tmp_dir}/input/{field} {output_dir}/{filter}-{field}-{input}.tsv  "
                    f" --codec zstd lz4 gzip blosclz --level 1 2 3 4 5 6 7 8 9 --filter {filter} --replicate {rep}"
                f" && rm -r {tmp_dir}",
                file=task_txt
            )

# Jobs for BPCells
for rep in replicates:
    output_dir = f"{RESULTS_ROOT}/raw/compression/in-memory-compression/bpcells_{rep}/"
    os.makedirs(output_dir, exist_ok=True)
    for input in input_paths.keys():
        tmp_dir = f"{TMP_ROOT}/in-memory-compression_{tmp_dir_id}"
        tmp_dir_id += 1
        print(
            f"mkdir -p {tmp_dir} "
            f" && {SINGULARITY} Rscript decompress-bpcells.R {input_paths[input]} {input_type[input]} {tmp_dir}/input"
            f" && {SINGULARITY} Rscript compress_{input_type[input]}_BPCells.R {tmp_dir}/input {rep} {output_dir}/{input}.tsv"
            f" && rm -r {tmp_dir}",
            file=task_txt
        )