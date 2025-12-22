import anndata as ad

adata = ad.read_h5ad("plate2_LungCancer_CVCL_0459.h5ad", backed='r')

# 将 obs 保存成 CSV（只有细胞元信息，内存可控）
adata.obs.to_csv("tahoe100M_plate2_LungCancer_CVCL_0459_obs.csv")
