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

# We have 23 subset sets with 10 minute time limits each (can go up to 2x time limit)
SBATCH_OPTS=" ".join([
    "-p biochem",
    "--constraint=CPU_SKU:7502",
    f"-c {THREADS}",
    "--time=5:00:00",
    "--mem=128G"
])

print(SBATCH_OPTS, file=task_txt)

CENSUS_ROOT = f"{DATA_ROOT}/cellxgene-census/2024-07-01/"

UNIQUE_CELLS = int(
    subprocess.run(
        ['Rscript', '-e', f'cat(length(readRDS("{CENSUS_ROOT}/obs_idx.rds")))'],
        capture_output=True
    ).stdout)
GENES = int(
    subprocess.run(
        ['Rscript', '-e', f'cat(nrow(readRDS("{CENSUS_ROOT}/var.rds")))'],
        capture_output=True
    ).stdout)

# Create slice sets (axis_name, axis_len, subset_sizes)
subset_info = [
    ("cell", UNIQUE_CELLS, [1, 10, 100, 1_000, 10_000, 100_000, 300_000]),
    ("gene", GENES, [1, 10, 100, 300]),
]

SLICE_VARIANTS=10
for axis, axis_len, subset_sizes in subset_info:
    for subset_size in subset_sizes:
        subset_dir = f"{CENSUS_ROOT}/slice_coords/{axis}_{subset_size}"
        if not os.path.exists(f"{subset_dir}/sequential_{SLICE_VARIANTS:03d}.txt"):
            print(f"Generating subset coords in dir: {subset_dir}")
            os.makedirs(subset_dir, exist_ok=True)
            subprocess.run(f"{SINGULARITY} Rscript {file_dir}/generate_slices.R {axis_len} {subset_size} {SLICE_VARIANTS} {subset_dir}", shell=True)

# Add an all-item subset
axis, axis_len, subset_size = ("gene", GENES, GENES)
subset_dir = f"{CENSUS_ROOT}/slice_coords/{axis}_{subset_size}"
if not os.path.exists(f"{subset_dir}/sequential_001.txt"):
    print(f"Generating subset coords in dir: {subset_dir}")
    os.makedirs(subset_dir, exist_ok=True)
    subprocess.run(f"{SINGULARITY} Rscript {file_dir}/generate_slices.R {axis_len} {subset_size} 1 {subset_dir}", shell=True)

if IS_LAPTOP:
    replicates = ["rep1", "rep2"]
else:
    replicates = ["rep1", "rep2", "rep3"]

# Input format, slice command
slice_targets = [
    ("bpcells_cellmajor", "Rscript slice_bpcells.R"),
    ("bpcells_genemajor", "Rscript slice_bpcells.R"),
    ("bpcells_cellmajor_norm", "Rscript slice_bpcells.R"),
    ("bpcells_genemajor_norm", "Rscript slice_bpcells.R"),
    ("tiledb", "python slice_tiledb.py"),
    ("tiledb_norm", "python slice_tiledb.py"),
]

# Steps
# - For each subset size, run time_limited_runner.py with 10 minutes time limit
#    - This will loop over the candidate subset selections, breaking the loop after >10 minutes have elapsed
MAX_MINUTES = 10

os.makedirs(f"{RESULTS_ROOT}/raw/cellxgene/02_matrix_slicing/", exist_ok=True)
for (tool, command), rep in itertools.product(slice_targets, replicates):
    # Prep inputs once for each tool/replicate combination, then run all subsets on the same node
    staging_dir = f"{TMP_ROOT}/matrix_slicing_{tool}_{rep}"
    
    commands = []
    commands.append(f"mkdir -p {staging_dir}")
    commands.append(f"{SINGULARITY} Rscript prep_inputs.R {tool} {CENSUS_ROOT} {staging_dir} {IS_LAPTOP}")

    for axis, _, subset_sizes in subset_info:
        # Only run a full dataset read for the normalized variants
        if "_norm" in tool:
            if axis == "cell":
                continue
            else:
                subset_sizes = []
        if axis == "gene":
            subset_sizes = subset_sizes + [GENES]
        for subset_size in subset_sizes:
            output_dir = f"{RESULTS_ROOT}/raw/cellxgene/02_matrix_slicing/{tool}_{rep}/{axis}_{subset_size}"
            os.makedirs(output_dir, exist_ok=True)
            subset_dir = f"{CENSUS_ROOT}/slice_coords/{axis}_{subset_size}"
            
            # Run timing on a particular subset (where {{subset}} will get swapped out for its name)
            command_template = f"env time -vo {output_dir}/{{subset}}.gnutime.txt {command} {staging_dir} {subset_dir}/{{subset}}.txt {axis} {output_dir}/{{subset}}.tsv {THREADS}"
            commands.append(f"{SINGULARITY} python time_limited_runner.py '{command_template}' {subset_dir} {MAX_MINUTES}")
            commands.append(f"echo Starting {axis}_{subset_size} && date")
    commands.append(f"rm -r {staging_dir}")
    print(" && ".join(commands), file=task_txt)



