import anndata as ad
import pandas as pd
import numpy as np

# ======= 配置路径 =======
train_path = "competition_support_set/competition_train.h5"   # 训练数据
test_csv = "../vcc_data/pert_counts_Test.csv"     # 比赛提供的测试 perturbation 列表
output_path = "competition_support_set/preprocessed_test_data.h5ad"  # 输出模板路径

# ======= 1️⃣ 载入训练数据 =======
print("Loading training data...")
adata = ad.read_h5ad(train_path)

# 检查 control 条件
control_label = "non-targeting"   # 通常 starter.toml 里写了 control_pert = "non-targeting"
if "target_gene" not in adata.obs:
    raise KeyError("❌ 未找到 obs['target_gene']，请确认列名是否正确。")

# 取出对照细胞池
control_cells = adata[adata.obs["target_gene"] == control_label]
print(f"Found {control_cells.n_obs} control cells.")

# ======= 2️⃣ 载入测试任务列表 =======
test_df = pd.read_csv(test_csv)
print(f"Loaded {len(test_df)} perturbations from test set.")

# ======= 3️⃣ 构建模板 =======
synthetic_cells = []
for i, row in test_df.iterrows():
    pert = row["target_gene"]
    n_cells = int(row["n_cells"])

    # 若测试要求的细胞数 > 对照池大小，则循环采样（允许重复）
    sampled = control_cells[np.random.choice(control_cells.n_obs, size=n_cells, replace=True)].copy()
    sampled.obs["target_gene"] = pert  # 改标签为该 perturbation
    
    synthetic_cells.append(sampled)

# 合并所有伪造 perturbation 的细胞
synthetic_adata = ad.concat(synthetic_cells, axis=0)
synthetic_adata.obs.reset_index(drop=True, inplace=True)

# ======= 4️⃣ 生成 .h5ad 文件 =======
synthetic_adata.write_h5ad(output_path)
print(f"✅ Saved test input template to {output_path}")
print(f"Total cells: {synthetic_adata.n_obs}, genes: {synthetic_adata.n_vars}")
