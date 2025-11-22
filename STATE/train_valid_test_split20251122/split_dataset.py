'''
包含non-targeting 控制组的三个数据集划分代码
注意文件路径
'''
import numpy as np
import anndata as ad

print("\n读取 anndata (.h5ad)...")
adata = ad.read_h5ad("../preprocessed_data/preprocessed_training_data_2000_with_obs.h5ad")
print(f"adata: {adata}")
# 提取所有 target genes（去除 non-targeting）
all_targets = adata.obs['target_gene'].unique().tolist()
all_targets = [g for g in all_targets if g != "non-targeting"]
ntc_adata = adata[adata.obs["target_gene"] == "non-targeting"] # 添加 non-targeting 控制组

print("扰动基因数量：", len(all_targets))  # 应为 150

# 固定随机种子以保证可复现
np.random.seed(42)

# 打乱顺序
np.random.shuffle(all_targets)

# 划分比例，可根据需要调整
val_ratio = 0.3
test_ratio = 0.2

n = len(all_targets)
n_val = int(n * val_ratio)
n_test = int(n * test_ratio)

val_genes = all_targets[:n_val]
test_genes = all_targets[n_val:n_val + n_test]
train_genes = all_targets[n_val + n_test:]

print("Validation genes:", len(val_genes))
print("Test genes:", len(test_genes))
print("Train genes:",len(train_genes))
# -------------------------
# 生成 validation AnnData
# -------------------------
adata_val = adata[adata.obs['target_gene'].isin(val_genes)].copy()
adata_test = adata[adata.obs['target_gene'].isin(test_genes)].copy()
adata_train=adata[adata.obs['target_gene'].isin(train_genes)].copy()

# Append the non-targeting controls to the example anndata if they're missing
adata_val=ad.concat([adata_val, ntc_adata])
adata_test=ad.concat([adata_test, ntc_adata])
adata_train=ad.concat([adata_train, ntc_adata])

# 保存文件
date="1122"
adata_val.write_h5ad(f"validation_set_{date}.h5ad")
adata_test.write_h5ad(f"test_set_{date}.h5ad")
adata_train.write_h5ad(f"training_set_{date}.h5ad")

print(f"保存完成：validation_set_{date}.h5ad, test_set_{date}.h5ad, training_set_{date}.h5ad")
