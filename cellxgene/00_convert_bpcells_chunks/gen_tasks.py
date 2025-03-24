import itertools
import math
from pathlib import Path
import random


file_dir = Path(__file__).parent.resolve()
task_txt = open(file_dir / "tasks.txt", "w")

root_dir = Path(__file__).parent.resolve()
while not (root_dir / "config_vars.sh").is_file() and not root_dir.parent == root_dir:
    root_dir = root_dir.parent
config = {}
exec(open(root_dir / "config_vars.sh").read(), config)

THREADS = 3
SINGULARITY = config["CONTAINER_TEMPLATE"].format(
    CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS
)
DATA_ROOT = config["DATA_ROOT"]

INPUT_TILEDB=f"{DATA_ROOT}/cellxgene-census/2024-07-01/homo_sapiens/"
OUT_PATH=f"{DATA_ROOT}/cellxgene-census/2024-07-01"

# Taken from https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html
total_cells = 74_322_510

task_txt = open(Path(__file__).parent.resolve() / "tasks.txt", "w")

# Batches of 10 should be fine finishing within 30 minutes, but leave some margin
# convert_metadata takes >32G of RAM apparently
SBATCH_OPTS = " ".join(
    [
        "-p wjg,sfgf,biochem,owners",
        "--mem=64G",
        f"-c {THREADS}",
        "--time=40:00",
    ]
)
print(SBATCH_OPTS, file=task_txt)

APPLY_BATCHES = True
BATCH_SIZE = 10


# Matrix conversion + transpose jobs
tasks = []
for i in range(math.ceil(total_cells / 100_000)):
    tasks.append(
        f"{SINGULARITY} Rscript convert_bpcells_chunks.R {INPUT_TILEDB} {OUT_PATH}/bpcells {i+1}"
    )


# Copy demo of itertools batched implementation for python <3.12
def batched(iterable, n):
    # batched('ABCDEFG', 3) → ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(itertools.islice(it, n)):
        yield batch


if APPLY_BATCHES:
    random.Random(151243).shuffle(tasks)
    for task_group in batched(tasks, n=BATCH_SIZE):
        print(" && ".join(task_group), file=task_txt)
else:
    for task in tasks:
        print(task, file=task_txt)

# Metadata conversion
print(
    f"{SINGULARITY} Rscript convert_metadata.R {INPUT_TILEDB} {DATA_ROOT}/reference_metadata {OUT_PATH}",
    file=task_txt,
)
