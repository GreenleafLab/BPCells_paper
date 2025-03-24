import sys
import scanpy as sc

assert len(sys.argv) == 3
input = sys.argv[1]
output = sys.argv[2]

adata = sc.read_h5ad(input)
adata.var_names = adata.var.GeneID.values
adata = adata[adata.obs.method=="3GEX",:]
adata.write(output)