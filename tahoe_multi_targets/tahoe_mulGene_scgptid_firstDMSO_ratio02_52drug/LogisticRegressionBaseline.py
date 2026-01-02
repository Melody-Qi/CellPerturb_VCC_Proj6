#!/usr/bin/env python
"""
Logistic Regression baseline:
Input: delta expression per (drug, cell_line) + cell_line one-hot
Output: per-gene target probability (0/1 label)
"""

import ast
import json
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from tqdm import tqdm

# -----------------------------
# Config
# -----------------------------
H5AD_PATH = "h5ad/tahoe_log1p.h5ad"
DRUG_SUMMARY_CSV = "drug_summary.csv"
GENE_META_CSV = "gene_metadata.csv"
DRUG_SPLIT_JSON = "tahoe_drug_split_seed42.json"
SEED = 42
DEBUG = True

def debug(*args):
    if DEBUG:
        print("[DEBUG]", *args)

# -----------------------------
# Load metadata
# -----------------------------
debug("Loading gene metadata...")
gene_meta = pd.read_csv(GENE_META_CSV)
symbol2token = dict(zip(gene_meta["gene_symbol"], gene_meta["token_id"]))

debug("Loading drug summary...")
drug_df = pd.read_csv(DRUG_SUMMARY_CSV)

def parse_list(x):
    if pd.isna(x):
        return []
    if isinstance(x, list):
        return x
    try:
        return list(ast.literal_eval(x))
    except:
        return []

drug_df["target_genes"] = drug_df["target_genes"].apply(parse_list)
drug_df["cell_line_ids"] = drug_df["cell_line_ids"].apply(parse_list)

# -----------------------------
# Load drug split
# -----------------------------
with open(DRUG_SPLIT_JSON) as f:
    split_json = json.load(f)
drug_split = split_json["drug_split"]
train_drugs = drug_split["train"]
val_drugs = drug_split["val"]
test_drugs = drug_split["test"]

# -----------------------------
# Open h5ad
# -----------------------------
debug("Opening h5ad...")
adata = sc.read_h5ad(H5AD_PATH, backed="r")
var_tokens = adata.var.index.astype(int).values
token2idx = {tok:i for i, tok in enumerate(var_tokens)}
debug("adata shape:", adata.shape)

# -----------------------------
# Get all cell lines for one-hot encoding
# -----------------------------
all_cell_lines = sorted({cl for ids in drug_df["cell_line_ids"] for cl in ids})
cl2idx = {cl:i for i, cl in enumerate(all_cell_lines)}
n_cell_lines = len(all_cell_lines)
debug("n_cell_lines:", n_cell_lines)

# -----------------------------
# Build X, y for a list of drugs
# -----------------------------
def build_Xy(drug_list):
    X_list = []
    y_list = []
    info_list = []

    for drug in tqdm(drug_list, desc="Building X/y"):
        row = drug_df[drug_df["drug"]==drug].iloc[0]
        target_tokens = [symbol2token[g] for g in row["target_genes"] if g in symbol2token]
        if len(target_tokens)==0:
            debug(drug,"skipped: no target genes in gene_meta")
            continue

        for cl in row["cell_line_ids"]:
            # select cells
            mask = (adata.obs["drug"]==drug) & (adata.obs["cell_line_id"]==cl)
            if mask.sum()==0:
                debug(drug, cl, "skipped: no cells in h5ad")
                continue

            sub = adata[mask].to_memory()
            delta = sub.X.mean(axis=0) - sub.layers["ctrl_norm"].mean(axis=0)
            delta = np.asarray(delta).ravel()

            # cell_line one-hot
            cl_vec = np.zeros(n_cell_lines)
            cl_vec[cl2idx[cl]] = 1.0

            # concat delta + cell_line one-hot
            X_vec = np.concatenate([delta, cl_vec])
            X_list.append(X_vec)

            # y: 1 for target genes, 0 else
            y_vec = np.zeros_like(delta)
            for tok in target_tokens:
                idx = token2idx.get(tok)
                if idx is not None:
                    y_vec[idx] = 1
            y_list.append(y_vec)
            info_list.append({"drug":drug,"cell_line":cl})

    if len(X_list)==0:
        return None, None, None

    X = np.vstack(X_list)          # shape: n_samples x (n_genes + n_cell_lines)
    y = np.vstack(y_list)          # shape: n_samples x n_genes
    return X, y, info_list

# -----------------------------
# Build train / val / test
# -----------------------------
X_train, y_train, info_train = build_Xy(train_drugs)
X_val, y_val, info_val = build_Xy(val_drugs)
X_test, y_test, info_test = build_Xy(test_drugs)

# -----------------------------
# Flatten for logistic regression
# -----------------------------
# Each row = (drug, cell_line)
# Each column = gene delta+one-hot
n_samples, n_features = X_train.shape
n_genes = y_train.shape[1]
debug("n_samples:", n_samples, "n_features:", n_features, "n_genes:", n_genes)

X_train_flat = np.repeat(X_train, n_genes, axis=0)  # repeat per gene
y_train_flat = y_train.flatten()
X_val_flat   = np.repeat(X_val, n_genes, axis=0)
y_val_flat   = y_val.flatten()
X_test_flat  = np.repeat(X_test, n_genes, axis=0)
y_test_flat  = y_test.flatten()

# Add gene index as additional feature
gene_idx_feature = np.tile(np.arange(n_genes), X_train.shape[0]).reshape(-1,1)
X_train_flat = np.hstack([X_train_flat, gene_idx_feature])
X_val_flat   = np.hstack([X_val_flat, gene_idx_feature])
X_test_flat  = np.hstack([X_test_flat, gene_idx_feature])

# -----------------------------
# Train logistic regression
# -----------------------------
clf = LogisticRegression(max_iter=200, random_state=SEED, solver="liblinear")
debug("Training logistic regression...")
clf.fit(X_train_flat, y_train_flat)

# -----------------------------
# Evaluate
# -----------------------------
def evaluate(X, y, info_list, set_name="val"):
    y_pred_prob = clf.predict_proba(X)[:,1]
    roc = roc_auc_score(y, y_pred_prob)
    pr = average_precision_score(y, y_pred_prob)
    print(f"{set_name} ROC-AUC: {roc:.4f}, PR-AUC: {pr:.4f}")
    # optional: save predictions per drug/cell_line
    pred_df = pd.DataFrame({
        "drug":[info_list[i//n_genes]["drug"] for i in range(len(y))],
        "cell_line":[info_list[i//n_genes]["cell_line"] for i in range(len(y))],
        "gene_idx":np.tile(np.arange(n_genes), len(info_list)),
        "y_true":y,
        "y_pred_prob":y_pred_prob
    })
    pred_df.to_csv(f"logreg_{set_name}_pred.csv", index=False)
    return pred_df

evaluate(X_val_flat, y_val_flat, info_val, "val")
evaluate(X_test_flat, y_test_flat, info_test, "test")
