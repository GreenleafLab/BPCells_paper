import os
from pathlib import Path

file_dir = Path(__file__).parent.resolve()
task_txt = open(file_dir / "tasks.txt", "w")

root_dir = Path(__file__).parent.resolve()
while not (root_dir / "config_vars.sh").is_file() and not root_dir.parent == root_dir:
    root_dir = root_dir.parent
config = {}
exec(open(root_dir / "config_vars.sh").read(), config)

# Note: limiting time factor is import into ArchR / Snapatac2. 
# ArchR seems to require >16 GB RAM for 500k_heart sample ENCSR051ECW
# ArchR requires about 80 minutes for 1m_brain sample OCG_2 
# Fragment merge of 500k_heart normally runs in ~30 minutes, but weirdly in some cases
#   it appears to run extremely slowly (>3 hours). I'm not sure if this is a bad node or something
THREADS=4
SBATCH_OPTS=" ".join([
    "-p wjg,biochem,sfgf,owners",
    f"-c {THREADS}",
    "--time=3:00:00",
    "--mem=24G"
])
print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-locked-v1", OMP_NUM_THREADS=1)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"


# Consolidate the cell-barcodes files for datasets where we had to manually calculate per-sample
for dataset in ["35k_hematopoiesis", "45k_pancreas", "500k_heart"]:
    if not Path(f"{DATA_ROOT}/atac/{dataset}/cell-barcodes.txt").exists():
        print(f"Consolidating the per-sample barcode whitelist for {dataset}")
        cell_barcodes = []
        for sample in Path(f"{DATA_ROOT}/atac/{dataset}/cell-barcodes/").iterdir():
            assert sample.suffixes == [".txt"]
            cell_barcodes += open(sample).readlines()
        open(f"{DATA_ROOT}/atac/{dataset}/cell-barcodes.txt", "w").writelines(cell_barcodes)

# Merge and filter fragment files for each dataset
datasets = ["1m_brain", "3k_pbmc", "10k_pbmc", "35k_hematopoiesis", "45k_pancreas", "500k_heart"]
for dataset in datasets:
    print(f"{SINGULARITY} Rscript atac-filter-merge.R {DATA_ROOT}/atac/{dataset} {THREADS}", file=task_txt)

# Calculate peak and tile matrices for each dataset
for dataset in datasets:
    print(f"{SINGULARITY} Rscript calculate-matrices.R {DATA_ROOT}/atac/{dataset} {THREADS} {REFERENCE_DIR}", file=task_txt)

# Convert to ArchR arrow files
for dataset in datasets:
    os.makedirs(f"{DATA_ROOT}/atac/{dataset}/archr", exist_ok=True)
    for fragment_file in Path(f"{DATA_ROOT}/atac/{dataset}/fragments").iterdir():
        if ".fragments.tsv.gz.tbi" in fragment_file.name:
            continue
        sample = fragment_file.name.replace(".fragments.tsv.gz", "")
        print(f"{SINGULARITY} Rscript import-archr.R {DATA_ROOT}/atac/{dataset} {sample}", file=task_txt)

# Convert to snapatac2 h5ad files (separate calls for filtered and unfiltered)
for dataset in datasets:
    os.makedirs(f"{DATA_ROOT}/atac/{dataset}/snapatac2", exist_ok=True)
    for fragment_file in Path(f"{DATA_ROOT}/atac/{dataset}/fragments").iterdir():
        if ".fragments.tsv.gz.tbi" in fragment_file.name:
            continue
        sample = fragment_file.name.replace(".fragments.tsv.gz", "")
        print(f"{SINGULARITY} python import-snapatac2.py {DATA_ROOT}/atac/{dataset} {sample} {REFERENCE_DIR}", file=task_txt)
        print(f"{SINGULARITY} python import-snapatac2.py {DATA_ROOT}/atac/{dataset} {sample} {REFERENCE_DIR} --filter", file=task_txt)