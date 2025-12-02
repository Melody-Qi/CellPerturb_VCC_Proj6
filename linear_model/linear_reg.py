#!/usr/bin/env python3
"""
linear_reg.py

Usage example:
python linear_reg.py --input-h5ad ../vcc_data/small_set.h5ad --pred-out ./pred/small_pred_lr.h5ad --model-out ./pred/lr_baseline.joblib --ridge-alpha 1.0
python linear_reg.py --input-h5ad ../vcc_data/training_set_1119.h5ad --pred-out ./pred/train_pred_lr.h5ad --model-out ./pred/lr_baseline.joblib --ridge-alpha 1.0
"""
import argparse
import os
import numpy as np
import anndata as ad
import scanpy as sc
from scipy import sparse
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import Ridge
import joblib

def safe_normalize_log1p(adata):
    # 在副本上做 normalize_total + log1p，返回新 AnnData（不覆盖输入）
    a = ad.AnnData(X=adata.X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    if sparse.issparse(a.X):
        a.X = a.X.astype(np.float32)
    else:
        a.X = np.asarray(a.X, dtype=np.float32)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a

def build_design_matrix(adata, cat_cols=["target_gene","guide_id","batch"]):
    obs_df = adata.obs[cat_cols].astype(str).copy()
    enc = OneHotEncoder(sparse_output=True, handle_unknown='ignore')
    X_sparse = enc.fit_transform(obs_df.values)  # sparse (n_cells, n_features)
    feature_names = []
    for i, col in enumerate(cat_cols):
        cats = enc.categories_[i]
        feature_names += [f"{col}__{c}" for c in cats]
    print("Design matrix shape:", X_sparse.shape, "nnz:", X_sparse.nnz)
    for i, cats in enumerate(enc.categories_):
        print(f"  {cat_cols[i]} categories: {len(cats)}")
    return X_sparse.tocsr(), enc, np.array(feature_names)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--pred-out", required=True)
    parser.add_argument("--model-out", required=False)
    parser.add_argument("--ridge-alpha", type=float, default=1.0)
    args = parser.parse_args()

    print("Loading AnnData:", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    print("adata shape:", adata.shape, "sparse:", sparse.issparse(adata.X), "dtype:", adata.X.dtype)

    # 1) Preprocess: normalize + log1p (training space)
    print("Applying normalize_total + log1p to adata (creating a new adata_log1p)...")
    adata_log1p = safe_normalize_log1p(adata)

    # 2) Prepare Y (dense) -- WARNING: can be large. print memory hint.
    print("Converting adata_log1p.X to dense Y (float32). This may use significant memory.")
    if sparse.issparse(adata_log1p.X):
        Y = adata_log1p.X.toarray().astype(np.float32)
    else:
        Y = np.asarray(adata_log1p.X, dtype=np.float32)
    print("Y shape:", Y.shape, "dtype:", Y.dtype, "nbytes ~", Y.nbytes)

    # 3) Design matrix X (sparse) from categorical covariates
    X_sparse, enc, feature_names = build_design_matrix(adata_log1p, cat_cols=["target_gene","guide_id","batch"])

    # 4) Fit Ridge (multioutput)
    # Use solver that's okay for sparse design matrices; 'sparse_cg' works with sparse X
    print("Fitting Ridge(alpha={}) ...".format(args.ridge_alpha))
    model = Ridge(alpha=args.ridge_alpha, fit_intercept=True, solver='sparse_cg')
    # sklearn Ridge supports sparse X for fit when solver='sparse_cg' or 'lsqr'
    model.fit(X_sparse, Y)  # multi-output
    print("Fitted. coef_ shape:", model.coef_.shape, "intercept_ shape:", model.intercept_.shape)

    # 5) Predict on training set (same X)
    print("Predicting on training set...")
    Y_pred = model.predict(X_sparse)  # dense array
    Y_pred = Y_pred.astype(np.float32)

    # 6) Ensure no negative values (log1p space shouldn't be negative, but numerical eps may)
    Y_pred = np.maximum(Y_pred, 0.0)

    # 7) Build pred AnnData and save
    print("Constructing pred AnnData and saving to:", args.pred_out)
    pred_adata = ad.AnnData(X=Y_pred, obs=adata.obs.copy(), var=adata.var.copy())
    pred_dir = os.path.dirname(args.pred_out)
    if pred_dir and not os.path.exists(pred_dir):
        os.makedirs(pred_dir, exist_ok=True)
    pred_adata.write_h5ad(args.pred_out)
    print("Saved predicted AnnData to", args.pred_out)

    # 8) Save model + encoder + metadata
    if args.model_out:
        save_obj = {
            "model": model,
            "onehot_encoder": enc,
            "feature_names": feature_names,
            "var_names": adata.var_names.values,
            "preproc_use_log1p": True,
            "ridge_alpha": args.ridge_alpha
        }
        joblib.dump(save_obj, args.model_out)
        print("Saved model+encoder to", args.model_out)

    print("Done.")

    # adata = ad.read_h5ad("../vcc_data/validation_set_1119.h5ad")
    # valid_log1p=safe_normalize_log1p(adata)
    # X_sparse, enc, feature_names = build_design_matrix(valid_log1p, cat_cols=["target_gene","guide_id","batch"])
    # Y_pred = model.predict(X_sparse)  # dense array
    # Y_pred = Y_pred.astype(np.float32)
    # Y_pred = np.maximum(Y_pred, 0.0)
    # pred_adata = ad.AnnData(X=Y_pred, obs=adata.obs.copy(), var=adata.var.copy())
    # pred_dir = os.path.dirname("./pred/valid_pred_lr.h5ad")
    # if pred_dir and not os.path.exists(pred_dir):
    #     os.makedirs(pred_dir, exist_ok=True)
    # pred_adata.write_h5ad("./pred/valid_pred_lr.h5ad")
    # print("Saved predicted AnnData to", "./pred/valid_pred_lr.h5ad")

if __name__ == "__main__":
    main()
