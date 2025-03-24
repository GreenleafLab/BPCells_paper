import argparse
from pathlib import Path

import snapatac2

parser = argparse.ArgumentParser(description="Import sample to snapatac2 format")
parser.add_argument("dataset_dir", type=Path, help="Dataset folder")
parser.add_argument("sample_id", type=str, help="Sample ID")
parser.add_argument("reference_dir", type=Path, help="Path of folder with reference files")
parser.add_argument("--filter", action='store_true', help="If true, restrict to importing whitelist cells")
args = parser.parse_args()

genome = open(args.dataset_dir / "genome.txt").read().strip()
   
keeper_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
chrom_sizes_all = dict(tuple(l.strip().split('\t')) for l in open(args.reference_dir / (genome + ".chrom.sizes")))
chrom_sizes = {c: int(chrom_sizes_all[c]) for c in keeper_chromosomes}

if args.filter:
    out_path = args.dataset_dir / f"snapatac2/{args.sample_id}.filtered.h5ad"
    cell_barcodes = [l.strip() for l in open(args.dataset_dir / "cell-barcodes.txt")]
    if open(args.dataset_dir / "cell-barcodes-prefix.txt").read().strip() == "true":
        prefix = f"{args.sample_id}."
        cell_barcodes = [c.removeprefix(prefix) for c in cell_barcodes if c.startswith(prefix)]
    whitelist = cell_barcodes
else:
    out_path = args.dataset_dir / f"snapatac2/{args.sample_id}.h5ad"
    whitelist = None

adata = snapatac2.pp.import_data(
    args.dataset_dir / f"fragments/{args.sample_id}.fragments.tsv.gz",
    chrom_sizes,
    file = out_path,
    whitelist = whitelist,
    min_num_fragments=0,
    sorted_by_barcode=False,
    shift_right=1
)

