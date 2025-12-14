import scanpy as sc
import numpy as np
import scipy.sparse as sp
from tqdm import tqdm
import anndata as ad
import gc


def safe_normalize_log1p(file_path):
    # "test_set_1119.h5ad"
    a = ad.read_h5ad(f"raw/{file_path}.h5ad")
    print("Before:", type(a.X), a.X.dtype, a.X.max())
    #  normalize_total + log1p
    if sp.issparse(a.X):
        a.X = a.X.astype(np.float32)
    else:
        a.X = np.asarray(a.X, dtype=np.float32)
    # sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    print("After:", type(a.X), a.X.dtype, a.X.max())
    a.write_h5ad(f"log1p/{file_path}_lognorm.h5ad")
    print("Saved to log1p/..._lognorm.h5ad")


input_path = "../vcc_data/raw/test.h5ad"   # ←← 修改为你的原始数据路径
output_path = "../vcc_data/raw/test_baseline.h5ad" # baseline 输出文件（英文名）
non_target_name = "non-targeting"

print("Loading input h5ad...")
adata = sc.read_h5ad(input_path)
X = adata.X.tocsr().astype(np.float32)
targets = adata.obs["target_gene"].values
n_obs, n_genes = X.shape

print(f"Shape: {n_obs} x {n_genes}, nnz={X.nnz}")

#1) compute mean per target (exclude non-targeting)
print("Computing mean expression per target_gene (excluding non-targeting)...")
unique_targets = np.unique(targets)
target_mean = {}

for tg in tqdm(unique_targets):
    if tg == non_target_name:
        continue
    idx = np.where(targets == tg)[0]
    mean_row = X[idx].mean(axis=0)
    mean_row = np.asarray(mean_row).ravel().astype(np.float32)
    target_mean[tg] = sp.csr_matrix(mean_row)

#2) keep only non-targeting original rows, then free full X
nt_idx = np.where(targets == non_target_name)[0]
print(f"Copying non-targeting rows: {len(nt_idx)} cells ...")
X_nt = X[nt_idx].tocsr()  # only keep these original rows

# Free the big original X
del X
gc.collect()

# 3) Two-pass CSR build with a pointer into X_nt
print("Pass-1: counting nnz per row...")
indptr = np.zeros(n_obs + 1, dtype=np.int64)

j = 0  # pointer for non-targeting rows in X_nt
for i, tg in enumerate(tqdm(targets)):
    if tg == non_target_name:
        nnz = X_nt[j].nnz
        j += 1
    else:
        nnz = target_mean[tg].nnz
    indptr[i + 1] = indptr[i] + nnz

total_nnz = int(indptr[-1])
print(f"Total nnz for baseline: {total_nnz}")

print("Allocating CSR arrays...")
data = np.empty(total_nnz, dtype=np.float32)
indices = np.empty(total_nnz, dtype=np.int32)

print("Pass-2: filling CSR arrays...")
cursor = 0
j = 0
for i, tg in enumerate(tqdm(targets)):
    if tg == non_target_name:
        row = X_nt[j]
        j += 1
    else:
        row = target_mean[tg]

    nnz = row.nnz
    data[cursor:cursor + nnz] = row.data
    indices[cursor:cursor + nnz] = row.indices
    cursor += nnz

X_baseline = sp.csr_matrix((data, indices, indptr), shape=(n_obs, n_genes), dtype=np.float32)

print("Writing baseline.h5ad ...")
adata_out = sc.AnnData(X=X_baseline, obs=adata.obs.copy(), var=adata.var.copy())
adata_out.uns["baseline_type"] = "target_mean_baseline_non_targeting_identity"
adata_out.write(output_path)

print("Done:", output_path)