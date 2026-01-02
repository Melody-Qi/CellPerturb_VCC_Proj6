#!/usr/bin/env python
"""
Test whether drug target genes appear in top-K delta genes
using parquet_cache/global_parquet_stats.pkl
"""

import pickle
import ast
import numpy as np
import pandas as pd
from tqdm import tqdm

# -----------------------------
# Config
# -----------------------------
GLOBAL_STATS_PKL = "parquet_cache/global_parquet_stats.pkl"
DRUG_SUMMARY_CSV = "drug_summary.csv"
GENE_META_CSV = "gene_metadata.csv"

TOPKS = [5, 10, 50, 100]
DEBUG = True

def debug(*args):
    if DEBUG:
        print("[DEBUG]", *args)

# -----------------------------
# Load gene metadata
# -----------------------------
print("Loading gene metadata...")
gene_meta = pd.read_csv(GENE_META_CSV)
symbol2token = dict(zip(gene_meta["gene_symbol"], gene_meta["token_id"]))
debug(f"Loaded {len(symbol2token)} gene symbols")

# -----------------------------
# Load drug summary
# -----------------------------
print("Loading drug summary...")
drug_df = pd.read_csv(DRUG_SUMMARY_CSV)

def parse_gene_list(x):
    if pd.isna(x):
        return []
    if isinstance(x, list):
        return x
    try:
        return list(ast.literal_eval(x))
    except Exception:
        return []

drug_df["target_genes"] = drug_df["target_genes"].apply(parse_gene_list)

# drug -> target genes
drug2targets = dict(
    zip(drug_df["drug"], drug_df["target_genes"])
)

# -----------------------------
# Load global parquet stats
# -----------------------------
print("Loading global parquet stats...")
with open(GLOBAL_STATS_PKL, "rb") as f:
    global_stats = pickle.load(f)

print(f"Loaded {len(global_stats)} (cell_line, drug) entries")

# -----------------------------
# Core logic
# -----------------------------
results = []

for (cell_line, drug), stats in tqdm(global_stats.items(), desc="Processing (cell_line, drug)"):
    delta = stats["delta"]            # dict: token -> delta
    total_count = stats["total_count"]

    target_genes = drug2targets.get(drug, [])
    if not target_genes:
        continue

    # map target genes to tokens
    target_tokens = {}
    for g in target_genes:
        tok = symbol2token.get(g)
        if tok is not None:
            target_tokens[g] = tok

    if not target_tokens:
        continue

    # -------------------------
    # rank genes by |delta|
    # -------------------------
    tokens = np.array(list(delta.keys()))
    deltas = np.array([delta[t] for t in tokens])

    order = np.argsort(np.abs(deltas))[::-1]
    ranked_tokens = tokens[order]

    token2rank = {tok: i + 1 for i, tok in enumerate(ranked_tokens)}

    # -------------------------
    # per-target statistics
    # -------------------------
    gene_ranks = {}
    topk_counts = {k: 0 for k in TOPKS}

    for g, tok in target_tokens.items():
        rank = token2rank.get(tok, np.inf)
        gene_ranks[g] = rank
        for k in TOPKS:
            if rank <= k:
                topk_counts[k] += 1

    record = {
        "drug": drug,
        "cell_line_id": cell_line,
        "n_rows": total_count,
        "target_genes": ",".join(target_genes),
        "gene_ranks": gene_ranks,
        "best_rank": min(gene_ranks.values())
    }

    for k in TOPKS:
        record[f"n_in_top{k}"] = topk_counts[k]

    results.append(record)

# -----------------------------
# Save results
# -----------------------------
res_df = pd.DataFrame(results)
res_df.to_csv("drug_target_in_topk_from_parquet_cache.csv", index=False)

print("\n✅ Done!")
print("Total records:", len(res_df))
print(res_df.head())
