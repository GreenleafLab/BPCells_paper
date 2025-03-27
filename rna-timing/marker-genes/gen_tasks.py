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


# Expand our node list to include all SH3_CBASE configs across the whole cluster,
# since they should be built from identical hardware
THREADS=16
SBATCH_OPTS=" ".join([
    "-p biochem,owners",
    "--constraint=CPU_SKU:7502",
    "--constraint=CLASS:SH3_CBASE",
    f"-c {THREADS}",
    "--time=3:00:00",
    "--mem=128000M"
])

print(SBATCH_OPTS, file=task_txt)

# Make a second task.txt for the big-memory jobs
os.makedirs(file_dir / "high_mem_tasks", exist_ok = True)
task_txt_highmem = open(file_dir / "high_mem_tasks" / "tasks.txt", "w")
print(SBATCH_OPTS.replace("128000M", "256000M"), file=task_txt_highmem)

DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"
SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)

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
small_datasets = ["20k_pbmc", "130k_thymus_atlas", "480k_tabula_sapiens", "500k_drugscreen", "1m_neurons"]

# These datasets are excluded because tools can't finish within time + memory constraints
dataset_exclusions = {
    "bpcells": [],
    # R dgCMatrix with Seurat can't handle > 2^31 non-zero entries. We'll let it go up to 1m to prove the crash
    "presto": ["2m_perturbseq", "4m_fetal", "11m_jax", "22m_pansci"],
    # We should almost certainly get timeouts or memory crashes before these 
    "scanpy": ["11m_jax", "22m_pansci"],
    # I expect timeouts for everything 130k cells and later
    "scanpy_tiecorrect": ["1m_neurons", "2m_perturbseq", "4m_fetal", "11m_jax", "22m_pansci"],
}

# These datasets should be run using the maximum available memory
dataset_high_memory = {
    "bpcells": [],
    "presto": [],
    # TODO: check if Scanpy really needs all this memory 
    "scanpy": ["1m_neurons", "2m_perturbseq", "4m_fetal", "11m_jax", "22m_pansci"],
    "scanpy_tiecorrect": [],
}

replicates = ["rep1","rep2","rep3","rep4", "rep5"]
tools = ["presto", "scanpy", "scanpy_tiecorrect", "bpcells"]

commands = {
    "bpcells": "Rscript markers_bpcells.R",
    "presto": "Rscript markers_presto.R",
    "scanpy": "python markers_scanpy.py notcorrected",
    "scanpy_tiecorrect": "python markers_scanpy.py tiecorrect",
}

# Group together all the BPCells thread counts for the small datasets, since they should run quickly
# and we can avoid spawning a bunch of jobs for the different thread counts
for rep, dataset, tool in itertools.product(replicates, datasets, tools):
    dataset_dir = f"{DATA_ROOT}/rna/{dataset}"
    output_dir = f"{RESULTS_ROOT}/raw/rna-timing/marker-genes/{tool}_{dataset}/"
    output_prefix = f"{output_dir}/{rep}"
    data_dir = f"{DATA_ROOT}/rna/marker-genes/{tool}"
    tmp_dir = f"{TMP_ROOT}/marker-genes_{tool}_{dataset}_{rep}"
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    command = (
        f"mkdir -p {tmp_dir}"
        f" && {SINGULARITY} Rscript prep_inputs.R {dataset_dir} {tmp_dir} {tool}"
        f" && {SINGULARITY} env time -vo {output_prefix}.gnutime.txt {commands[tool]} {tmp_dir} {output_prefix}.tsv {data_dir}/{dataset}_{rep}.csv"
        f" && rm -r {tmp_dir}"
    )

    if dataset in dataset_exclusions[tool]:
        command = f"echo 'mkdir -p {tmp_dir}' && echo Skipping dataset due to not running within time/memory constraints"

    if dataset in dataset_high_memory[tool]:
        print("cd .. && " + command, file=task_txt_highmem)
    else:
        print(command, file=task_txt)
