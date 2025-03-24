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
# These datasets can be grouped together when BPCells is running, as they don't take much time to analyze
small_datasets = ["20k_pbmc", "130k_thymus_atlas", "480k_tabula_sapiens", "500k_drugscreen", "1m_neurons"]

# These datasets should be run using the maximum available memory
dataset_high_memory = {
    # DelayedArray times out on anything big enough to cause memory worries
    "delayedarray": [],
    # R dgCMatrix with Seurat can't handle > 2^31 non-zero entries
    "seurat": [],
    # Scanpy crashes OOM with 128GB allocation with 1m and up 
    # (though it should be able to do at least 1m with up to 256GB RAM)
    "scanpy": ["1m_neurons", "2m_perturbseq", "4m_fetal", "11m_jax", "22m_pansci"],
    # Scanpy_dask has OOM crashes consistently for the largest 3 datasets
    "scanpy_dask": ["4m_fetal", "11m_jax", "22m_pansci"],
    "bpcells_stage_int": [],
    "bpcells_stage_float": [],
    "bpcells_stage_none": [],
    "bpcells_randomized": [],
}

replicates = ["rep1","rep2","rep3","rep4", "rep5"]
tools = ["delayedarray", "seurat", "scanpy", "bpcells_stage_int", "bpcells_stage_float", "bpcells_stage_none", "scanpy_dask", "bpcells_randomized"]

commands = {
    "delayedarray": "Rscript run_delayedarray.R",
    "seurat": "Rscript run_seurat.R",
    "scanpy": "python run_scanpy.py",
    "scanpy_dask": "python run_scanpy_dask.py",
    "bpcells_stage_int": "Rscript run_bpcells.R stage_int",
    "bpcells_stage_float": "Rscript run_bpcells.R stage_float",
    "bpcells_stage_none": "Rscript run_bpcells.R stage_none",
    "bpcells_randomized": "Rscript run_bpcells.R randomized",
}

thread_counts = {
    "delayedarray": [1, 2],
    "seurat": [1, 2],
    "scanpy": [1], # Scanpy doesn't support multiple cores for these workflows
    "scanpy_dask": [1,2,4,8,16],
    "bpcells_stage_int": [1,2,4,8,16],
    "bpcells_stage_float": [1,2,4,8,16],
    "bpcells_stage_none": [1,2,4,8,16],
    "bpcells_randomized": [1,2,4,8,16],
}

# These specific (tool, sample, thread) combinations appear to need high memory
extra_high_memory = [
    ("seurat", "500k_drugscreen", 2),
    ("scanpy_dask", "1m_neurons", 16),
    ("scanpy_dask", "2m_perturbseq", 8),
    ("scanpy_dask", "2m_perturbseq", 16),
]
extra_high_memory_tasks = []

# Group together all the BPCells thread counts for the small datasets, since they should run quickly
# and we can avoid spawning a bunch of jobs for the different thread counts
for tool, sample in itertools.product(tools, datasets):
    for rep in replicates:
        bpcells_grouped_jobs = []
        for threads in thread_counts[tool]:
            SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=threads)

            results_dir = f"{RESULTS_ROOT}/raw/rna-timing/pca-benchmark/{tool}/{threads}__{sample}__{rep}"
            os.makedirs(results_dir, exist_ok=True)
            
            sample_dir = f"{DATA_ROOT}/rna/{sample}"
            data_dir = f"{DATA_ROOT}/rna/pca/{tool}/{threads}__{sample}__{rep}"
            os.makedirs(data_dir, exist_ok=True)

            tmp_dir = f"{TMP_ROOT}/pca_{tool}/{threads}__{sample}__{rep}"

            if tool == "scanpy_dask":
                scanpy_dask_copy = f" && {SINGULARITY} Rscript -e \"library(BPCells); open_matrix_dir('{sample_dir}/bpcells') |> convert_matrix_type('float') |> write_matrix_anndata_hdf5('{tmp_dir}/input.h5ad')\"" 
            else:
                scanpy_dask_copy = ""
            
            command = (
                f"mkdir -p {tmp_dir}"
                f"{scanpy_dask_copy}"
                f" && {SINGULARITY} env time -vo {results_dir}/normalize.gnutime.txt {commands[tool]} {sample_dir} {data_dir} {results_dir} {tmp_dir} normalize"
                f" && {SINGULARITY} env time -vo {results_dir}/pca.gnutime.txt       {commands[tool]} {sample_dir} {data_dir} {results_dir} {tmp_dir} pca"
                f" && rm -r {tmp_dir}"
            )

            if "bpcells" in tool and sample in small_datasets:
                bpcells_grouped_jobs.append(command)
            elif sample in dataset_high_memory[tool]:
                print("cd .. && " + command, file=task_txt_highmem)
            elif (tool, sample, threads) in extra_high_memory:
                extra_high_memory_tasks.append("cd .. && " + command)
            else:
                print(command, file=task_txt)
        if len(bpcells_grouped_jobs) > 0:
            print(" && ".join(bpcells_grouped_jobs),file=task_txt)

# Print the extra high memory tasks at the end to not disturb earlier task IDs
for task in extra_high_memory_tasks:
    print(task, file=task_txt_highmem)