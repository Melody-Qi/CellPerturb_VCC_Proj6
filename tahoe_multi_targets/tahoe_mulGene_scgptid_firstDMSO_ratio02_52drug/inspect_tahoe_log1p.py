import scanpy as sc

adata = sc.read_h5ad("h5ad/tahoe_log1p.h5ad")

print(adata)
