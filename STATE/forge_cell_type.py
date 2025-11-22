import anndata as ad

# 读取两个文件
adata_raw = ad.read_h5ad("./competition_support_set/competition_train.h5")
adata_hvg = ad.read_h5ad("./preprocessed_data/preprocessed_training_data_2000.h5ad")

# 把原始 obs 补回 HVG 数据
adata_hvg.obs["batch_var"] = adata_raw.obs["batch_var"]
adata_hvg.obs["cell_type"] = adata_raw.obs["cell_type"]

# 保存新的 h5ad
adata_hvg.write_h5ad("./preprocessed_data/preprocessed_training_data_2000_with_obs.h5ad")

print("Done! obs 现在包含 batch_var & cell_type")
