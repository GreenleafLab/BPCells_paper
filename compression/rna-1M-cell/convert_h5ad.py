import os.path
import sys
import time

import scanpy as sc

assert len(sys.argv) == 3
input_path = sys.argv[1]
output_path = sys.argv[2]

assert input_path.endswith(".h5")

m = sc.read_10x_h5(input_path)
m.write_h5ad(output_path)