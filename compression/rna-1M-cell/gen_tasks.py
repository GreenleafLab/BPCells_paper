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

# Try to use the whole machine, since these timings are pretty short and hence sensitive to
# interference on memory bandwidth, etc.
THREADS=30
SBATCH_OPTS=" ".join([
    "-p biochem,owners",
    "--constraint=CPU_SKU:7502",
    "--constraint=CLASS:SH3_CBASE",
    f"-c {THREADS}",
    "--time=0:20:00",
    "--mem=250G"
])


print(SBATCH_OPTS, file=task_txt)

SINGULARITY=config["CONTAINER_TEMPLATE"].format(CONTAINER="bpcells-v0.3.0", OMP_NUM_THREADS=THREADS)
DATA_ROOT=config["DATA_ROOT"]
TMP_ROOT=config["TMP_ROOT"]
RESULTS_ROOT=config["RESULTS_ROOT"]
REFERENCE_DIR = f"{DATA_ROOT}/reference_metadata"

replicates = ["rep1","rep2","rep3","rep4", "rep5"]

path_10x_in = f"{DATA_ROOT}/rna/downloads/1m_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5"

for rep in replicates:
    tmp_dir = f"{TMP_ROOT}/rna-1m_{rep}"
    output_dir = f"{RESULTS_ROOT}/raw/compression/rna-1M-cell/{rep}/"
    os.makedirs(output_dir, exist_ok=True)

    print(
        f"mkdir -p {tmp_dir}"
        # Setup input files
        f" && cp {path_10x_in} {tmp_dir}/10x.h5"
        f" && {SINGULARITY} python convert_h5ad.py {tmp_dir}/10x.h5 {tmp_dir}/anndata.h5ad"
        f" && {SINGULARITY} Rscript convert_bpcells.R {tmp_dir}/10x.h5 {tmp_dir}/bpcells"
        # Use sleep to try to get similar disk performance state between scanpy and bpcells runs
        f" && sleep 1m" 
        f" && {SINGULARITY} python time_scanpy.py {tmp_dir}/10x.h5 {tmp_dir} {output_dir}/scanpy.tsv "
        f" && sleep 1m " 
        f" && {SINGULARITY} Rscript time_bpcells.R {tmp_dir}/10x.h5 {tmp_dir}/anndata.h5ad {tmp_dir} {output_dir}/bpcells.tsv "
        # Disk copies to have some control data about disk speed
        f" && {SINGULARITY} env time -vo {output_dir}/cp_10x.gnutime.txt cp {tmp_dir}/10x.h5 {tmp_dir}/10x.cp "
        f" && {SINGULARITY} env time -vo {output_dir}/cp_h5ad.gnutime.txt cp {tmp_dir}/mat.h5ad {tmp_dir}/h5ad.cp "
        f" && rm -r {tmp_dir}",
        file=task_txt
    )
 