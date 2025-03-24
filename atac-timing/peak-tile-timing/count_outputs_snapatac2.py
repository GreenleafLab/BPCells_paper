import sys

import snapatac2


args = sys.argv[1:]
assert len(args) == 3

staging_dir = args[0]
matrix_type = args[1]
output_path = args[2]

assert matrix_type in ["peak", "tile"]

if matrix_type == "peak":
    adata = snapatac2.read(f"{staging_dir}/matrix.h5ad")
elif matrix_type == "tile":
    adata = snapatac2.read(f"{staging_dir}/sample.h5ad")


total_overlaps = sum(c[0].sum() for c in adata.X.chunked(1024))

open(output_path, "w").write(f"{total_overlaps}\n")