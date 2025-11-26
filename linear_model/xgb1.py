#!/usr/bin/env python3
"""
xgb1.py

Usage examples:
python xgb1.py --input-h5ad ../vcc_data/small_set.h5ad --pred_out ./pred/small_pred_xgb.h5ad --model_out ./pred/xgb_baseline.joblib --use-log1p
python xgb1.py --input-h5ad ../vcc_data/training_set_1119.h5ad --pred_out ./pred/training_xgb.h5ad --model_out ./pred/xgb_baseline.joblib --use-log1p

Outputs:
    - pred_out (h5ad) with predicted X, obs, var preserved.
    - optional model_out (joblib) containing model + encoder + metadata.
Notes:
    - Uses XGBoost via sklearn wrapper MultiOutputRegressor.
    - For large gene counts this trains one estimator per output internally; may be memory/time heavy.
"""
import argparse
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from sklearn.preprocessing import OneHotEncoder
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import r2_score, mean_absolute_error
import joblib
import xgboost as xgb
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------- helper functions ----------------
def load_and_preprocess(h5ad_path: str, use_log1p: bool = True, hvg=None):
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
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("Skipping log1p; using raw adata.X as Y")

    # optional HVG filtering (用户可以传入 list of gene names)
    if hvg is not None:
        print(f"Applying HVG filter: keeping {len(hvg)} genes")
        mask = np.isin(adata.var_names.values, np.array(hvg))
        adata = adata[:, mask].copy()
        print("adata shape after HVG:", adata.shape)

    return adata

def build_design_matrix(adata: ad.AnnData, cat_cols=["target_gene", "guide_id", "batch"]):
    """
    Build sparse one-hot design matrix from adata.obs categorical columns.
    Returns (X_sparse_csr, encoder, feature_names)
    """
    obs_df = adata.obs[cat_cols].astype(str).copy()
    enc = OneHotEncoder(sparse_output=True, handle_unknown='ignore')
    X_sparse = enc.fit_transform(obs_df.values)
    feature_names = []
    for i, col in enumerate(cat_cols):
        cats = enc.categories_[i]
        feature_names += [f"{col}__{c}" for c in cats]
    print("Design matrix shape:", X_sparse.shape, "; num features:", len(feature_names))
    return X_sparse.tocsr(), enc, np.array(feature_names)

def prepare_Y(adata: ad.AnnData):
    print("Preparing Y matrix from adata.X ... (this will convert to dense; watch memory)")
    if sparse.issparse(adata.X):
        Y = adata.X.toarray().astype(np.float32)
    else:
        Y = adata.X.astype(np.float32)
    var_names = adata.var['gene_id'].values if 'gene_id' in adata.var.columns else adata.var_names.values
    return Y, np.array(var_names)

def fit_xgboost_multioutput(X, Y, xgb_params=None, n_estimators=100, n_jobs=1, random_state=0, verbose=False):
    """
    Fit XGBoost regressor wrapped in MultiOutputRegressor.
    X: sparse csr or dense array, shape (n_samples, n_features)
    Y: dense array (n_samples, n_outputs)
    """
    print("Fitting XGBoost MultiOutputRegressor; X shape:", X.shape, "Y shape:", Y.shape)
    if xgb_params is None:
        xgb_params = {
            "objective": "reg:squarederror",
            "tree_method": "hist",
            "learning_rate": 0.1,
            "max_depth": 6,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "verbosity": 0,
            "random_state": random_state
        }
    # sklearn wrapper estimator for a single output
    base = xgb.XGBRegressor(n_estimators=n_estimators, **xgb_params)
    # MultiOutputRegressor will spawn one estimator per target; set n_jobs to parallelize
    mor = MultiOutputRegressor(base, n_jobs=n_jobs)
    # convert sparse to array if XGBoost wrapper cannot accept sparse in sklearn wrapper
    X_train = X.toarray() if sparse.issparse(X) else X
    mor.fit(X_train, Y)
    if verbose:
        print("Finished training XGBoost MultiOutputRegressor.")
    return mor

def evaluate_on_train(model, X, Y, var_names, top_k=5):
    # optional quick training diagnostics
    try:
        X_in = X.toarray() if sparse.issparse(X) else X
        Y_pred = model.predict(X_in)
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
    parser = argparse.ArgumentParser(description="Train XGBoost baseline and write prediction h5ad for evaluation.")
    parser.add_argument("--input-h5ad", required=True, help="Training AnnData h5ad (used to fit baseline).")
    parser.add_argument("--pred_out", required=True, help="Output prediction h5ad path (e.g. ../small_set_pred_xgb.h5ad).")
    parser.add_argument("--model_out", default=None, help="Optional path to save model+encoder (joblib).")
    parser.add_argument("--use-log1p", action="store_true", help="Normalize+log1p before training/prediction.")
    parser.add_argument("--no-eval", action="store_true", help="Skip small training diagnostics (default: do quick eval).")
    parser.add_argument("--n-estimators", type=int, default=100, help="Number of boosting rounds for each XGBRegressor.")
    parser.add_argument("--n-jobs", type=int, default=1, help="Number of parallel jobs for MultiOutputRegressor.")
    parser.add_argument("--learning-rate", type=float, default=0.1, help="XGBoost learning rate.")
    parser.add_argument("--max-depth", type=int, default=6, help="XGBoost max depth.")
    parser.add_argument("--hvg", default=None, help="Optional CSV file with gene list (one gene per line) to restrict training to HVGs.")
    parser.add_argument("--cat-cols", default="target_gene,guide_id,batch", help="Comma-separated categorical columns in obs to encode.")
    args = parser.parse_args()

    # load + preprocess
    hvg_list = None
    if args.hvg:
        if os.path.exists(args.hvg):
            hvg_list = pd.read_csv(args.hvg, header=None)[0].astype(str).tolist()
        else:
            print("Warning: provided hvg file not found, ignoring --hvg")
            hvg_list = None

    adata = load_and_preprocess(args.input_h5ad, use_log1p=args.use_log1p, hvg=hvg_list)
    Y, var_names = prepare_Y(adata)

    cat_cols = [c.strip() for c in args.cat_cols.split(",") if c.strip()]
    X, enc, feature_names = build_design_matrix(adata, cat_cols=cat_cols)

    # fit model
    xgb_params = {
        "objective": "reg:squarederror",
        "tree_method": "hist",
        "learning_rate": args.learning_rate,
        "max_depth": args.max_depth,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "verbosity": 0
    }
    model = fit_xgboost_multioutput(X, Y, xgb_params=xgb_params, n_estimators=args.n_estimators, n_jobs=args.n_jobs, random_state=42, verbose=True)

    if not args.no_eval:
        print("Running quick train evaluation (on training data)...")
        evaluate_on_train(model, X, Y, var_names, top_k=5)

    # Predict using trained model
    print("Predicting full expression matrix with trained model...")
    X_in = X.toarray() if sparse.issparse(X) else X
    Y_pred = model.predict(X_in)  # (n_cells, n_genes)
    Y_pred = Y_pred.astype(np.float32)
    # post-processing: negative clip, optional log1p is not applied here to match training flow
    Y_pred = np.maximum(Y_pred, 0)
    # cap extreme predictions to keep numeric stability (same as earlier script)
    Y_pred = np.minimum(Y_pred, 1e6)

    # Construct pred AnnData: preserve obs & var ordering
    print("Constructing pred AnnData and saving to:", args.pred_out if hasattr(args, "pred_out") else args.pred_out)
    pred_adata = ad.AnnData(X=Y_pred, obs=adata.obs.copy(), var=adata.var.copy())

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
            "preproc_use_log1p": args.use_log1p,
            "xgb_params": xgb_params,
            "n_estimators": args.n_estimators
        }
        joblib.dump(save_obj, args.model_out)
        print("Saved model + encoder to", args.model_out)

    print("Done. You can now run your evaluation script which expects true_file and pred_file.")

if __name__ == "__main__":
    main()
