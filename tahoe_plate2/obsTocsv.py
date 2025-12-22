import anndata as ad

adata = ad.read_h5ad("plate2_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad", backed='r')

# 将 obs 保存成 CSV（只有细胞元信息，内存可控）
adata.obs.to_csv("tahoe100M_obs.csv")
