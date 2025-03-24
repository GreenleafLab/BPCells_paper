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

# Note: We allocate >1 thread since snapatac2 uses multiple threads even if we don't ask for it
# This much memory is probably not required, but we'll allocate it to be safe (I expect <5GB for worst case)
THREADS=4
SBATCH_OPTS=" ".join([
    "-p biochem",
    "--constraint=CPU_SKU:7502",
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

datasets = ["1m_brain", "3k_pbmc", "10k_pbmc", "35k_hematopoiesis", "45k_pancreas", "500k_heart"]
replicates = ["rep1","rep2","rep3","rep4","rep5"]
peak_counts = [10, 1000, 100000]


# List of scripts to call for each step of benchmarking for each tool
peak_tools = {
    "bpcells": "Rscript peaks_bpcells.R",
    "archr": "Rscript peaks_archr.R",
    "snapatac2": "python peaks_snapatac2.py",
}
tile_tools = {
    "bpcells": "Rscript tiles_bpcells.R",
    "archr": "Rscript tiles_archr.R",
    "snapatac2": "python tiles_snapatac2.py"
}
copy_inputs_tools = {
    "bpcells": "bash copy_inputs_bpcells.sh",
    "archr": "bash copy_inputs_archr.sh",
    "snapatac2": "bash copy_inputs_snapatac2.sh",
}
count_outputs_tools = {
    "bpcells": "Rscript count_outputs_bpcells.R",
    "archr": "Rscript count_outputs_archr.R",
    "snapatac2": "python count_outputs_snapatac2.py",
}

# List of (dataset, sample) tuples
samples = []
for dataset in datasets:
    fragment_files = Path(f"{DATA_ROOT}/atac/{dataset}/fragments").glob("*.fragments.tsv.gz.tbi")
    samples.extend([(dataset, str(f.name).replace(".fragments.tsv.gz.tbi", "")) for f in fragment_files])

# Calculate peak subsets
for dataset in datasets:
    input_peaks = f"{DATA_ROOT}/atac/{dataset}/peaks.bed" 
    output_dir = f"{DATA_ROOT}/atac/{dataset}/peak-subsets/"
    chr_names_file = f"{DATA_ROOT}/atac/{dataset}/bpcells_filtered/merged/chr_names"
    os.makedirs(output_dir, exist_ok=True)
    for count in peak_counts:
        output_path = f"{output_dir}/{count}.bed"
        if os.path.exists(output_path):
            continue
        print(f"Generating {count}-peak subset for dataset {dataset}")
        subprocess.run(
            f"{SINGULARITY} Rscript subset_peaks.R {input_peaks} {output_path} {count} {chr_names_file}", shell=True
        )


# Calculate peak/tile matrices
# We have a standard set of scripts supplied for each tool that perform a set of standardized steps:
#  1. Copy data into the staging directory
#  2. Calculate the matrices
#  3. Calculate the sum of all values in the matrices
#  4. Delete the files from the temporary directory
# These scripts are split up so we can use the `time` command to track memory and CPU utilization
# of just the core matrix calculations

# Due to the extremely high number of tasks this generates (>12k), 
# group jobs 25 at a time which will cut us down to ~500 jobs likely running under an hour each (est < 40m worst case, 20m average).

APPLY_BATCHES=True
BATCH_SIZE=25

# Copy demo of itertools batched implementation for python <3.12
def batched(iterable, n):
    # batched('ABCDEFG', 3) → ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := tuple(itertools.islice(it, n)):
        yield batch

tmp_dir_counter = 0
for rep in replicates:
    tasks = []
    # Peak jobs
    for (dataset, sample), count, tool in itertools.product(samples, peak_counts, peak_tools):
        dataset_dir = f"{DATA_ROOT}/atac/{dataset}"
        output_dir = f"{RESULTS_ROOT}/raw/atac-timing/peak-tile-timing/{dataset}/peaks-{count}"
        output_prefix = f"{output_dir}/{sample}.{tool}.{rep}"
        peaks_file = f"{DATA_ROOT}/atac/{dataset}/peak-subsets/{count}.bed"
        tmp_dir = f"{TMP_ROOT}/{tmp_dir_counter}"
        tmp_dir_counter += 1
        os.makedirs(output_dir, exist_ok=True)
        tasks.append(
            f"{SINGULARITY} {copy_inputs_tools[tool]} {dataset_dir} {tmp_dir} {sample}"
            f" && {SINGULARITY} env time -vo {output_prefix}.gnutime.txt "
                    f"{peak_tools[tool]} {dataset_dir} {output_prefix}.tsv {peaks_file} {tmp_dir}"
            f" && {SINGULARITY} {count_outputs_tools[tool]} {tmp_dir} peak {output_prefix}.matrix_sum.txt"
            f" && rm -r {tmp_dir}"
        )
    # Tile jobs
    for (dataset, sample), tool in itertools.product(samples, peak_tools):
        dataset_dir = f"{DATA_ROOT}/atac/{dataset}"
        output_dir = f"{RESULTS_ROOT}/raw/atac-timing/peak-tile-timing/{dataset}/tiles"
        output_prefix = f"{output_dir}/{sample}.{tool}.{rep}"
        tmp_dir = f"{TMP_ROOT}/{tmp_dir_counter}"
        tmp_dir_counter += 1
        os.makedirs(output_dir, exist_ok=True)
        tasks.append(
            f"{SINGULARITY} {copy_inputs_tools[tool]} {dataset_dir} {tmp_dir} {sample}"
            f" && {SINGULARITY} env time -vo {output_prefix}.gnutime.txt "
                    f"{tile_tools[tool]} {dataset_dir} {output_prefix}.tsv {tmp_dir} {REFERENCE_DIR}"
            f" && {SINGULARITY} {count_outputs_tools[tool]} {tmp_dir} tile {output_prefix}.matrix_sum.txt"
            f" && rm -r {tmp_dir}"
        )
    if APPLY_BATCHES:
        random.Random(151243).shuffle(tasks)
        for task_group in batched(tasks, n=BATCH_SIZE):
            print(" && ".join(task_group),file=task_txt)
    else:
        for task in tasks:
            print(task, file=task_txt)