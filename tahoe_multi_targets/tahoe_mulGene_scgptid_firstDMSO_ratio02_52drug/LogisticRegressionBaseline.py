#!/usr/bin/env python
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import top_k_accuracy_score

# ----------------------------
# 配置
# ----------------------------
DATA_DIR = Path("data/tahoe_mulGene_52drug/h5ad")
K_LIST = [1, 5, 10, 50]  # Hit@K
RANDOM_SEED = 42

# ----------------------------
# 读取所有 h5ad 并合并
# ----------------------------
all_X = []
all_labels = []
all_groups = []  # cell_line_id 用于分组split

for file in DATA_DIR.glob("*.h5ad"):
    adata = sc.read_h5ad(file)
    X = adata.X
    if "ctrl" in adata.layers:
        ctrl = adata.layers["ctrl"]
    else:
        ctrl = np.zeros_like(X)
    
    # Δexpression
    delta_X = X - ctrl
    all_X.append(delta_X)
    
    # 标签
    labels = adata.obs["target_gene"].astype(str).values
    all_labels.extend(labels)
    
    # 分组
    all_groups.extend(adata.obs["cell_line_id"].astype(str).values)

    # 释放内存
    del adata

X = np.vstack(all_X)
y = np.array(all_labels)
groups = np.array(all_groups)

print(f"X shape: {X.shape}, y shape: {y.shape}")

# ----------------------------
# 标签编码
# ----------------------------
le = LabelEncoder()
y_encoded = le.fit_transform(y)
n_classes = len(le.classes_)
print(f"Number of target genes: {n_classes}")

# ----------------------------
# Train/Test Split (按 cell_line_id 分组)
# ----------------------------
gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=RANDOM_SEED)
train_idx, test_idx = next(gss.split(X, y_encoded, groups=groups))

X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

print(f"Train shape: {X_train.shape}, Test shape: {X_test.shape}")

# ----------------------------
# Logistic Regression (Multiclass)
# ----------------------------
clf = LogisticRegression(
    max_iter=1000,
    multi_class="multinomial",
    solver="lbfgs",
    n_jobs=-1,
    verbose=1,
    random_state=RANDOM_SEED
)

clf.fit(X_train, y_train)
y_prob = clf.predict_proba(X_test)

# ----------------------------
# 计算 Hit@K
# ----------------------------
print("\nHit@K Scores:")
for k in K_LIST:
    hit_k = top_k_accuracy_score(y_test, y_prob, k=k, labels=np.arange(n_classes))
    print(f"  Hit@{k}: {hit_k:.4f}")
