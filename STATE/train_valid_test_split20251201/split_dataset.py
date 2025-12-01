import pandas as pd
import anndata as ad
import numpy as np
import os
from sklearn.model_selection import train_test_split

# ================
# 设置随机种子
# ================
SEED = 42
np.random.seed(SEED)

# ================
# 1. 读入文件
# ================
TEST_GENES_PATH = "../../vcc_data/test_set.csv"
ADATA_PATH = "../preprocessed_data/preprocessed_training_data_2000_with_obs.h5ad"

df_test = pd.read_csv(TEST_GENES_PATH)
test_genes = set(df_test["Group"].astype(str).tolist())   # 助教给的30基因

adata = ad.read_h5ad(ADATA_PATH)

# perturbation 信息列
pert_col = "target_gene"    # 修改成你真实 obs 里的列名

all_genes = set(adata.obs[pert_col].unique())
all_genes = {g for g in all_genes if g != "non-targeting"}   # 去掉 control

print("总 perturbation 数：", len(all_genes))
print("助教 Test 集：", len(test_genes), "个基因")

# =========================
# 2. 获得 Train+Val 的可用基因
# =========================
trainval_genes = sorted(list(all_genes - test_genes))
print("Train+Val 可用基因数：", len(trainval_genes))

# =========================
# 3. 从120个基因中划出 Valid
# =========================
val_size = 0.2  # 20% → 大约24个基因
train_genes, val_genes = train_test_split(
    trainval_genes,
    test_size=val_size,
    random_state=SEED,
)

train_genes = set(train_genes)
val_genes = set(val_genes)

print("Train 基因：", len(train_genes))
print("Valid 基因：", len(val_genes))
print("Test 基因：", len(test_genes))

# =========================
# 4. 按基因划分细胞
# =========================
is_control = adata.obs[pert_col] == "non-targeting"

train_mask = adata.obs[pert_col].isin(train_genes) | is_control
val_mask   = adata.obs[pert_col].isin(val_genes)
test_mask  = adata.obs[pert_col].isin(test_genes)

adata_train = adata[train_mask].copy()
adata_val = adata[val_mask].copy()
adata_test = adata[test_mask].copy()

print("Train cells:", adata_train.n_obs)
print("Valid cells:", adata_val.n_obs)
print("Test cells:", adata_test.n_obs)

# =========================
# 5. 保存输出
# =========================
outdir = "."
os.makedirs(outdir, exist_ok=True)

adata_train.write_h5ad(f"{outdir}/train.h5ad")
adata_val.write_h5ad(f"{outdir}/val.h5ad")

# 保存 valid gene names
pd.DataFrame({"valid_genes": sorted(val_genes)}).to_csv(
    f"{outdir}/valid_genes.csv",
    index=False
)

# 保存 test gene names
with open(f"{outdir}/test_genes.txt", "w") as f:
    for g in sorted(test_genes):
        f.write(g + "\n")

print("\n✔ 数据划分完成，文件已生成:")
print(" - train.h5ad")
print(" - val.h5ad")
print(" - valid_genes.csv")
print(" - test_genes.txt")
print(" - 使用 SEED =", SEED)
