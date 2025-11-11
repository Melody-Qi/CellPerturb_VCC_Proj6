import scanpy as sc

adata = sc.read_h5ad("../vcc_data/preprocessed_training_data_2000.h5ad")
print(adata.obs.columns)
