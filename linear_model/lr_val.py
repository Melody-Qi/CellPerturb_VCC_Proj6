#!/usr/bin/env python3
"""
Evaluate lr_baseline.joblib on valid_set.h5ad and output prediction h5ad.
Produces the same format as training-time prediction h5ad.
"""

import numpy as np
import anndata as ad
import scanpy as sc
from scipy import sparse
import joblib
import os

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
    model_path = "pred/lr_baseline.joblib"
    valid_path = "../vcc_data/validation_set_1119.h5ad"
    pred_out = "pred/valid_pred_lr.h5ad"

    # load saved model
    print("Loading model:", model_path)
    saved = joblib.load(model_path)
    model = saved["model"]
    encoder = saved["onehot_encoder"]
    use_log1p = saved["preproc_use_log1p"]

    # load validation set
    adata = load_and_preprocess(valid_path, use_log1p=use_log1p)

    # prepare Y matrix (if you want metrics later) — not needed to produce h5ad
    if sparse.issparse(adata.X):
        Y = adata.X.toarray().astype(np.float32)
    else:
        Y = adata.X.astype(np.float32)

    # build design matrix
    X_valid = build_design_matrix_with_encoder(
        adata,
        encoder,
        cat_cols=["target_gene", "guide_id", "batch"]
    )

    # prediction
    print("Predicting...")
    Y_pred = model.predict(X_valid).astype(np.float32)

    # same post-processing as training script
    Y_pred = np.maximum(Y_pred, 0)
    Y_pred = np.log1p(Y_pred)

    # construct prediction AnnData with *same* obs/vars
    print("Building pred AnnData...")
    pred_adata = ad.AnnData(
        X=Y_pred,
        obs=adata.obs.copy(),
        var=adata.var.copy()
    )

    # save .h5ad
    out_dir = os.path.dirname(pred_out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    pred_adata.write_h5ad(pred_out)
    print("Saved prediction h5ad to:", pred_out)


if __name__ == "__main__":
    main()
