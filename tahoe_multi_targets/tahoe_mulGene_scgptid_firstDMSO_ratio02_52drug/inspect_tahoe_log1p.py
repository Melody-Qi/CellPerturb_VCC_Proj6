import scanpy as sc

adata = sc.read_h5ad("h5ad/tahoe_log1p.h5ad", backed="r")

print(adata)

print(adata.var.head())
