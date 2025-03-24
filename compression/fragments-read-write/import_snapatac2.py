from pathlib import Path
import os.path
import sys
import time

import snapatac2


args = sys.argv[1:]
assert len(args) == 4
dataset_dir = Path(args[0])
reference_dir = Path(args[1])
staging_dir = Path(args[2])
output_timing = args[3]

genome = open(dataset_dir / "genome.txt").read().strip()
keeper_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
chrom_sizes_all = dict(tuple(l.strip().split('\t')) for l in open(reference_dir / (genome + ".chrom.sizes")))
chrom_sizes = {c: int(chrom_sizes_all[c]) for c in keeper_chromosomes}

elapsed_start = time.time()
process_start = time.process_time()
adata = snapatac2.pp.import_data(
    staging_dir / "sample.fragments.tsv.gz",
    chrom_sizes,
    file = staging_dir / "sample.h5ad",
    min_num_fragments=0,
    sorted_by_barcode=False,
    shift_right = 1,
    n_jobs = 1
)
cpu_time = time.process_time() - process_start
elapsed_time = time.time() - elapsed_start

file_bytes = os.path.getsize(staging_dir / "sample.h5ad")

open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\tbytes\n",
    f"{cpu_time}\t{elapsed_time}\t{file_bytes}\n"
])