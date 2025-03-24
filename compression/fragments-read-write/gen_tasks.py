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
# Worst-case timing is expected to be over 1 hour for the 1m_brain OGC_2 sample,
# and with batched execution I expect to get up to 2.5h for a batch size of 10 
THREADS=4
SBATCH_OPTS=" ".join([
    "-p biochem",
    "--constraint=CPU_SKU:7502",
    f"-c {THREADS}",
    "--time=4:00:00",
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

# List of scripts to call for each step of benchmarking for each tool
import_tools = {
    "bpcells": "Rscript import_bpcells.R",
    "archr": "Rscript import_archr.R",
    "snapatac2": "python import_snapatac2.py",
    "10x": "python import_10x.py",
}
read_tools = {
    "bpcells": "Rscript read_bpcells.R",
    "archr": "Rscript read_archr.R",
    "snapatac2": "python read_snapatac2.py",
    "10x": "Rscript read_10x.R",
}

# List of (dataset, sample) tuples
samples = []
for dataset in datasets:
    fragment_files = Path(f"{DATA_ROOT}/atac/{dataset}/fragments").glob("*.fragments.tsv.gz.tbi")
    samples.extend([(dataset, str(f.name).replace(".fragments.tsv.gz.tbi", "")) for f in fragment_files])

# Calculate time + memory for import, followed by time to read from disk

# Due to the extremely high number of tasks this generates (>4k), 
# group jobs which will cut us down to ~500 jobs likely running (est. 3h worst-case)

APPLY_BATCHES=True
BATCH_SIZE=10

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
    for (dataset, sample), tool in itertools.product(samples, import_tools.keys()):
        dataset_dir = f"{DATA_ROOT}/atac/{dataset}"
        output_dir = f"{RESULTS_ROOT}/raw/compression/fragments-read-write/{dataset}/"
        output_prefix = f"{output_dir}/{sample}.{tool}.{rep}"
        tmp_dir = f"{TMP_ROOT}/{tmp_dir_counter}"
        tmp_dir_counter += 1
        os.makedirs(output_dir, exist_ok=True)
        tasks.append(
            f"{SINGULARITY} bash copy_inputs.sh {dataset_dir} {tmp_dir} {sample}"
            f" && {SINGULARITY} env time -vo {output_prefix}.import.gnutime.txt "
                    f"{import_tools[tool]} {dataset_dir} {REFERENCE_DIR} {tmp_dir} {output_prefix}.import.tsv "
            f" && {SINGULARITY} env time -vo {output_prefix}.read.gnutime.txt "
                    f"{read_tools[tool]} {dataset_dir} {REFERENCE_DIR} {tmp_dir} {output_prefix}.read.tsv "
            f" && rm -r {tmp_dir}"
        )
    
    if APPLY_BATCHES:
        random.Random(151243).shuffle(tasks)
        for task_group in batched(tasks, n=BATCH_SIZE):
            print(" && ".join(task_group),file=task_txt)
    else:
        for task in tasks:
            print(task, file=task_txt)