import anndata as ad
adata = ad.read_h5ad("submission/my_submission.vcc")
print(adata.uns.keys())
