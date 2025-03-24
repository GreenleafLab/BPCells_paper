
# We don't have to do real work here, just record file size in bytes
from pathlib import Path
import os.path
import sys
import time

args = sys.argv[1:]
assert len(args) == 4
dataset_dir = Path(args[0])
reference_dir = Path(args[1])
staging_dir = Path(args[2])
output_timing = args[3]

file_bytes = os.path.getsize(staging_dir / "sample.fragments.tsv.gz")

open(output_timing, "w").writelines([
    "time_cpu\ttime_elapsed\tbytes\n",
    f"NaN\tNaN\t{file_bytes}\n"
])