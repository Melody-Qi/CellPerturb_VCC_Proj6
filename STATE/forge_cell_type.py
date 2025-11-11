import anndata as ad

adata = ad.read_h5ad("../vcc_data/preprocessed_training_data_2000.h5ad")
adata.obs["cell_type"] = "H1"   # 全部设为同一个类型
adata.write_h5ad("../vcc_data/preprocessed_training_data_2000withcelltype.h5ad")
