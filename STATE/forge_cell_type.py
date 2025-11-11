import anndata as ad

import pandas as pd

adata = ad.read_h5ad("../vcc_data/preprocessed_training_data_2000.h5ad")
# adata.obs["cell_type"] = "H1"   # 全部设为同一个类型

adata.obs["cell_type"] = pd.Categorical(["H1"] * adata.n_obs)
adata.write_h5ad("../vcc_data/preprocessed_training_data_2000withcelltype.h5ad")
