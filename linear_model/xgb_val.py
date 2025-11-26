#!/usr/bin/env python3
"""
xgb_val.py

Usage:
python xgb_val.py --model-in pred/xgb_baseline.joblib --input-h5ad ../vcc_data/validation_set_1119.h5ad --pred-out pred/valid_pred_xgb.h5ad

Loads saved model+encoder and generate prediction h5ad.
"""

import argparse
import os
import numpy as np
import anndata as ad
import scanpy as sc
from scipy import sparse
import joblib

def load_and_preprocess(h5ad_path, use_log1p=True):
    print("Loading AnnData:", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    # ensure float32 CSR
    if sparse.issparse(adata.X):
        adata.X = adata.X.tocsr().astype(np.float32)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)

    if use_log1p:
        print("Normalizing + log1p...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    return adata

def build_design_matrix_with_encoder(adata, encoder, cat_cols=["target_gene","guide_id","batch"]):
    obs_df = adata.obs[cat_cols].astype(str).copy()
    X = encoder.transform(obs_df.values)  # DO NOT fit
    print("Design matrix shape:", X.shape)
    return X.tocsr()

def main():
    parser = argparse.ArgumentParser(description="Load XGBoost model joblib and predict on a new h5ad.")
    parser.add_argument("--model-in", required=True, help="Path to saved joblib model (contains model+encoder).")
    parser.add_argument("--input-h5ad", required=True, help="AnnData to predict on (e.g. validation set).")
    parser.add_argument("--pred-out", required=True, help="Output prediction h5ad path.")
    parser.add_argument("--cat-cols", default="target_gene,guide_id,batch", help="Comma-separated categorical columns in obs to encode.")
    args = parser.parse_args()

    print("Loading model:", args.model_in)
    saved = joblib.load(args.model_in)
    model = saved["model"]
    encoder = saved["onehot_encoder"]
    use_log1p = saved.get("preproc_use_log1p", True)

    adata = load_and_preprocess(args.input_h5ad, use_log1p=use_log1p)

    cat_cols = [c.strip() for c in args.cat_cols.split(",") if c.strip()]
    X_valid = build_design_matrix_with_encoder(adata, encoder, cat_cols=cat_cols)

    # prediction
    print("Predicting...")
    X_in = X_valid.toarray() if sparse.issparse(X_valid) else X_valid
    Y_pred = model.predict(X_in).astype(np.float32)

    # same post-processing as training script
    Y_pred = np.maximum(Y_pred, 0)
    Y_pred = np.minimum(Y_pred, 1e6)

    # construct prediction AnnData with *same* obs/vars
    print("Building pred AnnData...")
    pred_adata = ad.AnnData(
        X=Y_pred,
        obs=adata.obs.copy(),
        var=adata.var.copy()
    )

    out_dir = os.path.dirname(args.pred_out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    pred_adata.write_h5ad(args.pred_out)
    print("Saved prediction h5ad to:", args.pred_out)

if __name__ == "__main__":
    main()
