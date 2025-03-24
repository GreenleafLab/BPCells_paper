import os
from pathlib import Path
import urllib.request
import re
import subprocess
file_dir = Path(__file__).parent.resolve()
task_txt = open(file_dir / "tasks.txt", "w")

root_dir = Path(__file__).parent.resolve()
while not (root_dir / "config_vars.sh").is_file() and not root_dir.parent == root_dir:
    root_dir = root_dir.parent
config = {}
exec(open(root_dir / "config_vars.sh").read(), config)

# Note: limiting time factor is the 93 GB OGC_2 sample from catlas, running about 4 hours
# Without that, 3 hours is plenty of time
SBATCH_OPTS=" ".join([
    "-p wjg,biochem,sfgf,owners",
    "-c 4",
    "--time=6:00:00",
])
print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-locked-v1", OMP_NUM_THREADS=1)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]

# Download references (gencode, blacklists, chrom sizes) if they haven't already been downloaded
reference_urls = [
    'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz',
    'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz',
    'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz',
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"
reference_urls = [url for url in reference_urls 
                  if not os.path.exists(os.path.join(REFERENCE_DIR, os.path.basename(url)))]
os.makedirs(REFERENCE_DIR, exist_ok=True)
if len(reference_urls) > 0:
    print("Downloading reference files (e.g. Gencode). May take several minutes")
    args = ["curl", "--parallel"]
    for url in reference_urls:
        args.append(url)
        args.append("-o")
        args.append(os.path.join(REFERENCE_DIR, os.path.basename(url)))
    subprocess.run(args)


# 10x datasets
print(f"{SINGULARITY} bash 10x-datasets/download-10x-datasets.sh {DATA_ROOT}", file=task_txt)
print(f"{SINGULARITY} Rscript 10x-datasets/download_cell_ids.R {DATA_ROOT}", file=task_txt)


# catlas brain dataset
catlas_url = "http://catlas.org/catlas_downloads/humanbrain/bedpe/"
def get_catlas_samples(url_base):
    samples = []
    for l in urllib.request.urlopen(url_base).readlines():
        res = re.search(">([^<]+).bedpe.gz</a>", l.decode())
        # Note: CTXGA and PER have zero cells after filtering
        if res and res.group(1) not in ["CTXGA", "PER"]:
            samples.append(res.group(1))
    return samples
catlas_samples = get_catlas_samples(catlas_url)
dataset_dir = f"{DATA_ROOT}/atac/1m_brain"
for s in catlas_samples:
    fragment_path = f"{dataset_dir}/fragments/{s}.fragments.tsv.gz"
    bpcells_path = f"{dataset_dir}/bpcells/{s}"
    print(
        f"{SINGULARITY} bash brain/download_catlas_sample.sh {s} {TMP_ROOT} {fragment_path}",
        f"&& {SINGULARITY} Rscript fragments-to-bpcells.R {fragment_path} {bpcells_path}",
        file=task_txt
    )
print(f"{SINGULARITY} Rscript brain/get_cells.R {DATA_ROOT}/atac/1m_brain", file=task_txt)
catlas_peak_url="http://catlas.org/catlas_downloads/humanbrain/Supplementary_Tables/Table%20S6%20%e2%80%93%20List%20of%20cCREs%20in%20bed%20format.gz"
print(f"{SINGULARITY} bash -c 'curl \"{catlas_peak_url}\"| gunzip -c | cut -f 1-3 > {DATA_ROOT}/atac/1m_brain/peaks.bed'", file=task_txt)

# ENCODE heart and pancreas datasets
encode_heart_ids = [l.strip() for l in open(file_dir / "encode" / "ENCODE_scATAC_heart.txt")]
dataset_dir = f"{DATA_ROOT}/atac/500k_heart/"
os.makedirs(f"{dataset_dir}/cell-barcodes", exist_ok=True)
for s in encode_heart_ids:
    fragment_path = f"{dataset_dir}/fragments/{s}.fragments.tsv.gz"
    bpcells_path = f"{dataset_dir}/bpcells/{s}"
    passing_cells_path = f"{dataset_dir}/cell-barcodes/{s}.txt"
    print(
        f"{SINGULARITY} python encode/download_encode_sample.py {s} {TMP_ROOT} {fragment_path}", 
        f"&& {SINGULARITY} Rscript fragments-to-bpcells.R {fragment_path} {bpcells_path}",
        f"&& {SINGULARITY} Rscript find-passing-cells.R {bpcells_path} {passing_cells_path} hg38 4000 15 {REFERENCE_DIR}",
        file=task_txt
    )

encode_pancreas_ids = [l.strip() for l in open(file_dir / "encode" / "ENCODE_scATAC_pancreas.txt")]
dataset_dir = f"{DATA_ROOT}/atac/45k_pancreas/"
os.makedirs(f"{dataset_dir}/cell-barcodes", exist_ok=True)
for s in encode_pancreas_ids:
    fragment_path = f"{dataset_dir}/fragments/{s}.fragments.tsv.gz"
    bpcells_path = f"{dataset_dir}/bpcells/{s}"
    passing_cells_path = f"{dataset_dir}/cell-barcodes/{s}.txt"
    print(
        f"{SINGULARITY} python encode/download_encode_sample.py {s} {TMP_ROOT} {fragment_path}", 
        f"&& {SINGULARITY} Rscript fragments-to-bpcells.R {fragment_path} {bpcells_path}",
        f"&& {SINGULARITY} Rscript find-passing-cells.R {bpcells_path} {passing_cells_path} hg38 4000 15 {REFERENCE_DIR}",
        file=task_txt
    )

print(f"{SINGULARITY} bash encode/download_peakset.sh {DATA_ROOT}/atac", file=task_txt)

# 35k hematopoiesis dataset
print(f"{SINGULARITY} bash hematopoiesis/download-hematopoiesis.sh {TMP_ROOT} {DATA_ROOT}/atac/35k_hematopoiesis {REFERENCE_DIR}", file=task_txt)


# Write supplementary metadata to help later processing
dataset_metadata = {
    "3k_pbmc":          {"barcode_prefix": "false", "genome": "hg38"},
    "10k_pbmc":         {"barcode_prefix": "false", "genome": "hg38"},
    "35k_hematopoiesis":{"barcode_prefix": "true",  "genome": "hg19"},
    "45k_pancreas":     {"barcode_prefix": "true",  "genome": "hg38"},
    "500k_heart":       {"barcode_prefix": "true",  "genome": "hg38"},
    "1m_brain":         {"barcode_prefix": "false", "genome": "hg38"},
}
for d, v in dataset_metadata.items():
    dataset_dir=Path(f"{DATA_ROOT}/atac/{d}")
    os.makedirs(dataset_dir / "fragments", exist_ok=True)
    open(dataset_dir / "cell-barcodes-prefix.txt", "w").write(v["barcode_prefix"] + "\n")
    open(dataset_dir / "genome.txt", "w").write(v["genome"] + "\n")


