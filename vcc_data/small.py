import numpy as np
import anndata as ad

print("\n读取 anndata (.h5ad)...")
adata = ad.read_h5ad("adata_Training.h5ad")
all_targets = adata.obs['target_gene'].unique().tolist()
ntc_adata = adata[adata.obs["target_gene"] == "non-targeting"] 
small_genes = all_targets[:5]
adata_small = adata[adata.obs['target_gene'].isin(small_genes)].copy()
adata_small=ad.concat([adata_small, ntc_adata])
adata_small.write_h5ad(f"small_set.h5ad")

print(f"保存完成")