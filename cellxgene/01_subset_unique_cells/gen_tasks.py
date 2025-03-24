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


DATA_ROOT = config["DATA_ROOT"]
TMP_ROOT = config["TMP_ROOT"]
RESULTS_ROOT = config["RESULTS_ROOT"]
IS_LAPTOP = bool(config["IS_LAPTOP"])

# It appears that 5 hours is probably enough time for these jobs
if IS_LAPTOP:
    THREADS = 16
else:
    THREADS = 32
SBATCH_OPTS = " ".join(
    [
        "-p biochem",
        "--constraint=CPU_SKU:7502",
        f"-c {THREADS}",
        "--time=6:00:00",
        "--mem=250G",
    ]
)
SINGULARITY = config["CONTAINER_TEMPLATE"].format(
    CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS
)

print(SBATCH_OPTS, file=task_txt)

CENSUS_ROOT = f"{DATA_ROOT}/cellxgene-census/2024-07-01/"

if IS_LAPTOP:
    replicates = ["rep1", "rep2"]
else:
    replicates = ["rep1", "rep2", "rep3"]

# Handle tiledb subsetting
os.makedirs(f"{RESULTS_ROOT}/raw/cellxgene/01_subset_unique_cells/", exist_ok=True)
os.makedirs(f"{CENSUS_ROOT}/tiledb/", exist_ok=True)
for rep, compression_level, layer in itertools.product(
    replicates, ["1", "9"], ["normalized", "raw"]
):
    if IS_LAPTOP and compression_level != "9":
        continue
    tmp_dir = f"{TMP_ROOT}/tiledb_{layer}_zstd_{compression_level}_{rep}"
    input_path = f"{CENSUS_ROOT}/homo_sapiens/ms/RNA/X/raw"  # Always use raw input even if we're normalizing
    output_dir = f"{RESULTS_ROOT}/raw/cellxgene/01_subset_unique_cells/"
    output_prefix = f"{output_dir}/tiledb_{layer}_{compression_level}_{rep}"
    output_matrix = f"{CENSUS_ROOT}/tiledb/{layer}_zstd_{compression_level}"

    if IS_LAPTOP:
        buf_size_gb = 2
        base_command = f"mkdir -p {tmp_dir} && "
        command_input = f"{input_path}"
    else:
        buf_size_gb = 4
        base_command = f"mkdir -p {tmp_dir} && cp -r {input_path} {tmp_dir}/input && "
        command_input = f"{tmp_dir}/input"

    base_command = base_command + (
        f"{SINGULARITY} env time -vo {output_prefix}.gnutime.txt "
        f" python tiledb_write_{layer}.py {command_input} {tmp_dir}/output {CENSUS_ROOT} {compression_level} {output_prefix}.tsv {buf_size_gb}"
    )
    # Only save the outputs for replicate 1
    if rep == "rep1":
        command = base_command + f" && mv {tmp_dir}/output {output_matrix}"
    else:
        command = base_command
    command = command + f" && rm -r {tmp_dir}"

    print(command, file=task_txt)

UNIQUE_CELLS = int(
    subprocess.run(
        ["Rscript", "-e", f'cat(length(readRDS("{CENSUS_ROOT}/obs_idx.rds")))'],
        capture_output=True,
    ).stdout
)
OUTPUT_CHUNK_COUNT = math.ceil(UNIQUE_CELLS / 100_000)

# Handle BPCells subsetting
os.makedirs(f"{CENSUS_ROOT}/bpcells/", exist_ok=True)
for rep, orientation in itertools.product(
    replicates, ["cellmajor", "genemajor", "norm"]
):
    tmp_dir = f"{TMP_ROOT}/bpcells_{orientation}_{rep}"
    input_path = f"{CENSUS_ROOT}/bpcells/cellmajor_chunks"
    output_dir = f"{RESULTS_ROOT}/raw/cellxgene/01_subset_unique_cells/"
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = f"{output_dir}/bpcells_{orientation}_{rep}"
    output_matrix = f"{CENSUS_ROOT}/bpcells/subset_{orientation}"

    if IS_LAPTOP:
        base_command = f"mkdir -p {tmp_dir} && "
        command_input = f"{input_path}"
    else:
        base_command = f"mkdir -p {tmp_dir} && cp -r {input_path} {tmp_dir}/input && "
        command_input = f"{tmp_dir}/input"
    base_command = base_command + (
        f"{SINGULARITY} env time -vo {output_prefix}.gnutime.txt "
        f" Rscript bpcells_write_{orientation}.R {input_path} {CENSUS_ROOT} {OUTPUT_CHUNK_COUNT} {tmp_dir}/output {output_prefix}.tsv"
    )
    # Only save the outputs for replicate 1
    if rep == "rep1":
        command = base_command + f" && mv {tmp_dir}/output {output_matrix}"
    else:
        command = base_command
    command = command + f" && rm -r {tmp_dir}"

    print(command, file=task_txt)
