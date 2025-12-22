import anndata as ad
import numpy as np

# 1. 以 backed 模式打开（只读，不进内存）
adata = ad.read_h5ad(
    "plate2_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad",
    backed="r"
)

# 2. 选择 CVCL_0459
mask = adata.obs["cell_line"] == "CVCL_0459"

print("Selected cells:", np.sum(mask))

adata_sub = adata[mask].copy(
    filename="plate2_LungCancer_CVCL_0459.h5ad"
)
