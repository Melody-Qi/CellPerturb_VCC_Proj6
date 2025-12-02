#!/usr/bin/env python3
"""
eval_validation.py

Usage:
python lr_val.py --input-h5ad ../vcc_data/validation_set_1119.h5ad --model-in ./pred/lr_baseline.joblib --pred-out ./pred/valid_pred_lr.h5ad
"""

import argparse
import os
import numpy as np
import anndata as ad
import scanpy as sc
from scipy import sparse
import joblib


def safe_normalize_log1p(adata):
    """Same normalization as training."""
    a = ad.AnnData(X=adata.X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    if sparse.issparse(a.X):
        a.X = a.X.astype(np.float32)
    else:
        a.X = np.asarray(a.X, dtype=np.float32)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def build_design_matrix_with_existing_encoder(adata, encoder, cat_cols):
    """Use pre-fitted OneHotEncoder from training."""
    obs_df = adata.obs[cat_cols].astype(str).copy()
    X_sparse = encoder.transform(obs_df.values)
    return X_sparse.tocsr()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--model-in", required=True)
    parser.add_argument("--pred-out", required=True)
    args = parser.parse_args()

    print("Loading validation AnnData:", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    print("Shape:", adata.shape)

    print("Loading trained model + encoder:", args.model_in)
    saved = joblib.load(args.model_in)
    model = saved["model"]
    encoder = saved["onehot_encoder"]
    feature_names = saved["feature_names"]
    var_names_train = saved["var_names"]
    cat_cols = ["target_gene", "guide_id", "batch"]

    # 1) Normalize + log1p
    print("Applying normalization (normalize_total + log1p)...")
    adata_log1p = safe_normalize_log1p(adata)

    # 2) Build design matrix using training-time encoder
    print("Building design matrix...")
    X_sparse = build_design_matrix_with_existing_encoder(
        adata_log1p, encoder, cat_cols
    )
    print("Design matrix shape:", X_sparse.shape)

    # 3) Predict
    print("Predicting...")
    Y_pred = model.predict(X_sparse).astype(np.float32)
    Y_pred = np.maximum(Y_pred, 0.0)

    # 4) Construct output AnnData
    print("Constructing predicted AnnData and saving:", args.pred_out)

    # Make sure var order matches training model
    var_df = adata.var.copy()
    if not np.array_equal(var_df.index.values, var_names_train):
        print("WARNING: var_names mismatch; reindexing to training var order.")
        var_df = var_df.reindex(var_names_train)

    pred_adata = ad.AnnData(
        X=Y_pred,
        obs=adata.obs.copy(),
        var=var_df,
    )

    pred_dir = os.path.dirname(args.pred_out)
    if pred_dir and not os.path.exists(pred_dir):
        os.makedirs(pred_dir, exist_ok=True)
    pred_adata.write_h5ad(args.pred_out)

    print("Saved predicted AnnData to:", args.pred_out)
    print("Done.")


if __name__ == "__main__":
    main()
