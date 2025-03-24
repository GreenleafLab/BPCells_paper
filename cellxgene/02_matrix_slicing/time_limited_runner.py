
# Run a slice command template on each of the different subsets,
# but cut off after a maximum time limit 
# The command template is set as a string with '{subset}' marking where to substitute in the subset name

from pathlib import Path
import os.path
import time
import subprocess
import sys

args = sys.argv[1:]
assert len(args) == 3
subset_command_template = args[0]
coords_dir = Path(args[1])
maximum_minutes = int(args[2])

# Run timing tests 
sequential_subsets = sorted(coords_dir.glob(f"sequential_*.txt"))
random_subsets = sorted(coords_dir.glob(f"random_*.txt"))

end_time = time.time() + maximum_minutes * 60
for s in sequential_subsets:
    subprocess.run(subset_command_template.format(subset=s.stem), shell=True, check=True)
    if time.time() >= end_time:
        break

end_time = time.time() + maximum_minutes * 60
for s in random_subsets:
    subprocess.run(subset_command_template.format(subset=s.stem), shell=True, check=True)
    if time.time() >= end_time:
        break
