#!/usr/bin/env python
"""
Logistic Regression baseline using parquet_cache/global_parquet_stats.pkl

Each sample = (drug, cell_line, gene)
Feature = delta + cell_line one-hot
Label = whether gene is a target of the drug
"""

import pickle
import ast
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from tqdm import tqdm

# -----------------------------
# Config
# -----------------------------
GLOBAL_STATS_PKL = "parquet_cache/global_parquet_stats.pkl"
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

drug2targets = dict(zip(drug_df["drug"], drug_df["target_genes"]))
drug2celllines = dict(zip(drug_df["drug"], drug_df["cell_line_ids"]))

# -----------------------------
# Load drug split
# -----------------------------
with open(DRUG_SPLIT_JSON) as f:
    split_json = json.load(f)

train_drugs = split_json["drug_split"]["train"]
val_drugs   = split_json["drug_split"]["val"]
test_drugs  = split_json["drug_split"]["test"]

# -----------------------------
# Load global parquet stats
# -----------------------------
debug("Loading global parquet stats...")
with open(GLOBAL_STATS_PKL, "rb") as f:
    global_stats = pickle.load(f)

debug("Total (cell_line, drug):", len(global_stats))

# -----------------------------
# Collect all cell lines
# -----------------------------
all_cell_lines = sorted({
    cell for (cell, _) in global_stats.keys()
})
cl2idx = {cl: i for i, cl in enumerate(all_cell_lines)}
n_cell_lines = len(all_cell_lines)

debug("n_cell_lines:", n_cell_lines)

# -----------------------------
# Build dataset
# -----------------------------
def build_dataset(drug_list):
    X = []
    y = []

    for (cell_line, drug), stats in tqdm(global_stats.items()):
        if drug not in drug_list:
            continue

        target_genes = drug2targets.get(drug, [])
        if not target_genes:
            continue

        target_tokens = {
            symbol2token[g] for g in target_genes if g in symbol2token
        }
        if not target_tokens:
            continue

        delta = stats["delta"]  # dict: token -> delta

        # cell line one-hot
        cl_vec = np.zeros(n_cell_lines)
        cl_vec[cl2idx[cell_line]] = 1.0

        for tok, d in delta.items():
            # feature
            X_vec = np.concatenate([[d], cl_vec])
            X.append(X_vec)

            # label
            y.append(1 if tok in target_tokens else 0)

    X = np.vstack(X)
    y = np.array(y)
    return X, y

# -----------------------------
# Build splits
# -----------------------------
debug("Building train set...")
X_train, y_train = build_dataset(train_drugs)

debug("Building val set...")
X_val, y_val = build_dataset(val_drugs)

debug("Building test set...")
X_test, y_test = build_dataset(test_drugs)

debug("Train shape:", X_train.shape, "Positive rate:", y_train.mean())

# -----------------------------
# Train logistic regression
# -----------------------------
clf = LogisticRegression(
    max_iter=500,
    solver="liblinear",
    random_state=SEED
)

debug("Training logistic regression...")
clf.fit(X_train, y_train)

# -----------------------------
# Evaluate
# -----------------------------
def evaluate(X, y, name):
    prob = clf.predict_proba(X)[:, 1]
    roc = roc_auc_score(y, prob)
    pr  = average_precision_score(y, prob)
    print(f"{name} ROC-AUC: {roc:.4f}, PR-AUC: {pr:.4f}")

evaluate(X_val, y_val, "Val")
evaluate(X_test, y_test, "Test")

# =========================================================
# Per-drug gene ranking analysis on TEST set
# =========================================================

TOPKS = [5, 10, 50, 100]

print("\n=== Per-drug target gene ranking (TEST) ===")

# 重新跑一遍 test，用于收集 prediction + meta info
rows = []

for (cell_line, drug), stats in global_stats.items():
    if drug not in test_drugs:
        continue

    target_genes = drug2targets.get(drug, [])
    if not target_genes:
        continue

    target_tokens = {
        symbol2token[g] for g in target_genes if g in symbol2token
    }
    if not target_tokens:
        continue

    delta = stats["delta"]

    # cell line one-hot
    cl_vec = np.zeros(n_cell_lines)
    cl_vec[cl2idx[cell_line]] = 1.0

    for tok, d in delta.items():
        X_vec = np.concatenate([[d], cl_vec]).reshape(1, -1)
        prob = clf.predict_proba(X_vec)[0, 1]

        rows.append({
            "drug": drug,
            "cell_line": cell_line,
            "gene_token": tok,
            "prob": prob,
            "is_target": int(tok in target_tokens),
        })

pred_df = pd.DataFrame(rows)
print("Total test predictions:", len(pred_df))

# ---------------------------------------------------------
# Aggregate over cell lines: (drug, gene)
# ---------------------------------------------------------
agg_df = (
    pred_df
    .groupby(["drug", "gene_token"])
    .agg(
        mean_prob=("prob", "mean"),
        is_target=("is_target", "max")
    )
    .reset_index()
)

# ---------------------------------------------------------
# Rank genes per drug
# ---------------------------------------------------------
records = []

for drug, df_d in agg_df.groupby("drug"):
    df_d = df_d.sort_values("mean_prob", ascending=False).reset_index(drop=True)
    df_d["rank"] = np.arange(1, len(df_d) + 1)

    target_ranks = df_d[df_d["is_target"] == 1]["rank"].tolist()

    rec = {
        "drug": drug,
        "n_genes": len(df_d),
        "n_targets": len(target_ranks),
        "best_rank": min(target_ranks) if target_ranks else np.inf,
    }

    for k in TOPKS:
        rec[f"n_target_in_top{k}"] = sum(r <= k for r in target_ranks)

    records.append(rec)

    # 打印给你直观看
    print(f"\nDrug: {drug}")
    print("Target gene ranks:", target_ranks)
    for k in TOPKS:
        print(f"  Top-{k}: {sum(r <= k for r in target_ranks)}")

summary_df = pd.DataFrame(records)

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------
summary_df.to_csv("logreg_test_target_gene_ranking.csv", index=False)
agg_df.to_csv("logreg_test_gene_probs.csv", index=False)

print("\nSaved:")
print("  logreg_test_target_gene_ranking.csv")
print("  logreg_test_gene_probs.csv")

print("\n=== Summary ===")
print(summary_df)

