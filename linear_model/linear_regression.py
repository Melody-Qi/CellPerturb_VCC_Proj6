#!/usr/bin/env python3
"""
train_and_write_pred.py

Usage examples:
python linear_regression.py --input-h5ad ../vcc_data/small_set.h5ad --pred-out ./pred/small_pred_lr.h5ad --model-out ./pred/lr_baseline.joblib --use-log1p
python linear_regression.py --input-h5ad ../vcc_data/training_set_1119.h5ad --pred-out ./pred/training_lr.h5ad --model-out ./pred/lr_baseline.joblib --use-log1p

# memory-friendly (SGD incremental)
python linear_regression.py \
        --input-h5ad ../full_training.h5ad \
        --pred-out ../small_set_pred.h5ad \
        --use-sgd --chunk-size 20000

Outputs:
    - pred_out (h5ad) with predicted X, obs, var preserved.
    - optional model_out (joblib) containing model + encoder + metadata.
"""

import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import Ridge, SGDRegressor
from sklearn.metrics import r2_score, mean_absolute_error
import joblib
import os
import sys

# ---------------- helper functions ----------------
def load_and_preprocess(h5ad_path: str, use_log1p: bool = True):
    print("Loading AnnData:", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    print("adata shape:", adata.shape)
    # Ensure X is CSR float32
    if sparse.issparse(adata.X):
        adata.X = adata.X.tocsr().astype(np.float32)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)

    if use_log1p:
        print("Performing library-size normalization + log1p (inplace on adata.X)")
        # Work on adata directly; if user does not want, they can pass --no-log1p
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("Skipping log1p; using raw adata.X as Y")

    return adata

def build_design_matrix(adata: ad.AnnData, cat_cols=["target_gene", "guide_id", "batch"]):
    obs_df = adata.obs[cat_cols].astype(str).copy()
    enc = OneHotEncoder(sparse_output=True, handle_unknown='ignore')#把这些分类变量编码成稀疏矩阵（0/1）
    X_sparse = enc.fit_transform(obs_df.values) #(n_cells,n_features) 
    feature_names = []
    for i, col in enumerate(cat_cols):
        cats = enc.categories_[i]
        feature_names += [f"{col}__{c}" for c in cats]
    print("Design matrix shape:", X_sparse.shape, "; num features:", len(feature_names))
    return X_sparse.tocsr(), enc, np.array(feature_names)

def fit_ridge(X: sparse.csr_matrix, Y: np.ndarray, alpha: float = 1.0, random_state:int = 0):
    print("Fitting Ridge multi-output regression; X shape:", X.shape, "Y shape:", Y.shape)
    model = Ridge(alpha=alpha, fit_intercept=True, random_state=random_state)
    model.fit(X, Y)
    return model

def fit_sgd_incremental(X_sparse: sparse.csr_matrix, Y: np.ndarray, chunk_size: int = 20000,
                        alpha=1e-4, max_iter=1, random_state=0):
    """
    Memory-friendly incremental training using one SGDRegressor per gene.
    Warning: creates n_genes estimators; useful when dense solvers OOM.
    对于大数据集（例如百万级细胞），避免一次性加载全部数据。
    """
    n_samples, n_genes = Y.shape
    print(f"SGD incremental training: {n_samples} samples, {n_genes} genes, chunk_size={chunk_size}")
    estimators = [SGDRegressor(alpha=alpha, max_iter=max_iter, tol=None, warm_start=True, random_state=random_state+i)
                  for i in range(n_genes)]
    for start in range(0, n_samples, chunk_size):
        end = min(start + chunk_size, n_samples)
        print(f"  chunk {start}:{end}")
        X_chunk = X_sparse[start:end].toarray() if sparse.issparse(X_sparse) else X_sparse[start:end]
        Y_chunk = Y[start:end]
        for j, est in enumerate(estimators):
            est.partial_fit(X_chunk, Y_chunk[:, j])
    class Wrapper:
        def __init__(self, estimators, n_features):
            self.estimators = estimators
            self.n_features_in_ = n_features
        def predict(self, X):
            X_arr = X.toarray() if sparse.issparse(X) else X
            preds = np.column_stack([est.predict(X_arr) for est in self.estimators])
            return preds
    return Wrapper(estimators, X_sparse.shape[1])

def evaluate_on_train(model, X, Y, var_names, top_k=5):
    # optional quick training diagnostics
    try:
        Y_pred = model.predict(X)
    except Exception as e:
        print("Warning: model.predict failed during evaluation:", e)
        return None
    r2s = np.array([r2_score(Y[:, i], Y_pred[:, i]) for i in range(Y.shape[1])])
    maes = np.array([mean_absolute_error(Y[:, i], Y_pred[:, i]) for i in range(Y.shape[1])])
    df = pd.DataFrame({"gene": var_names, "r2": r2s, "mae": maes}).sort_values("r2", ascending=False)
    print("Train median R2:", np.nanmedian(r2s))
    print("Top genes by R2:\n", df.head(top_k))
    return df

# ---------------- main ----------------
def main():
    parser = argparse.ArgumentParser(description="Train linear baseline and write prediction h5ad for evaluation.")
    parser.add_argument("--input-h5ad", required=True, help="Training AnnData h5ad (used to fit baseline).")
    parser.add_argument("--pred-out", required=True, help="Output prediction h5ad path (e.g. ../small_set_pred.h5ad).")
    parser.add_argument("--model-out", default=None, help="Optional path to save model+encoder (joblib).")
    parser.add_argument("--use-log1p", action="store_true", help="Normalize+log1p before training/prediction.")
    parser.add_argument("--no-eval", action="store_true", help="Skip small training diagnostics (default: do quick eval).")
    parser.add_argument("--use-sgd", action="store_true", help="Use SGD incremental training (memory friendly).")
    parser.add_argument("--chunk-size", type=int, default=20000, help="Chunk size for SGD incremental.")
    parser.add_argument("--ridge-alpha", type=float, default=1.0, help="Alpha for Ridge.")
    args = parser.parse_args()

    adata = load_and_preprocess(args.input_h5ad, use_log1p=args.use_log1p)

    # Y matrix: convert to dense — NOTE: may be memory heavy for large datasets
    print("Preparing Y matrix from adata.X ... (this will convert to dense; watch memory)")
    if sparse.issparse(adata.X):
        Y = adata.X.toarray().astype(np.float32)
    else:
        Y = adata.X.astype(np.float32)
    var_names = adata.var['gene_id'].values if 'gene_id' in adata.var.columns else adata.var_names.values

    X, enc, feature_names = build_design_matrix(adata, cat_cols=["target_gene","guide_id","batch"])

    if args.use_sgd:
        model = fit_sgd_incremental(X, Y, chunk_size=args.chunk_size, alpha=1e-4, max_iter=1, random_state=42)
    else:
        model = fit_ridge(X, Y, alpha=args.ridge_alpha, random_state=42)

    if not args.no_eval:
        print("Running quick train evaluation (on training data)...")
        evaluate_on_train(model, X, Y, var_names, top_k=5)

    # Predict using trained model
    print("Predicting full expression matrix with trained model...")
    Y_pred = model.predict(X)  # (n_cells, n_genes)
    # ensure dtype float32
    Y_pred = Y_pred.astype(np.float32)
    Y_pred = np.maximum(Y_pred, 0)  # 负值截断为0
    #Y_pred = np.log1p(Y_pred)       # 再取log1p压缩
    #Y_pred = np.minimum(Y_pred, 1e5)

    # Construct pred AnnData: preserve obs & var ordering
    print("Constructing pred AnnData and saving to:", args.pred_out)
    pred_adata = ad.AnnData(X=Y_pred, obs=adata.obs.copy(), var=adata.var.copy())
    # keep same uns, obsm, varm if desired (here we drop by default)
    # Save
    pred_dir = os.path.dirname(args.pred_out)
    if pred_dir and not os.path.exists(pred_dir):
        os.makedirs(pred_dir, exist_ok=True)
    pred_adata.write_h5ad(args.pred_out)
    print("Saved predicted AnnData to", args.pred_out)

    # Optionally save model + encoder + metadata
    if args.model_out:
        save_obj = {
            "model": model,
            "onehot_encoder": enc,
            "feature_names": feature_names,
            "var_names": var_names,
            "preproc_use_log1p": args.use_log1p
        }
        joblib.dump(save_obj, args.model_out)
        print("Saved model + encoder to", args.model_out)

    print("Done. You can now run your evaluation script which expects true_file and pred_file.")

if __name__ == "__main__":
    main()
