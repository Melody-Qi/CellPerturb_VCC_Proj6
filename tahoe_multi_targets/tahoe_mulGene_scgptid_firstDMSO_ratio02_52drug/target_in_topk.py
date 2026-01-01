#!/usr/bin/env python
"""
Test whether drug target genes appear in top-K delta genes,
and report per-target gene ranking and top-K counts.
"""

import ast
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

# -----------------------------
# Config
# -----------------------------
H5AD_PATH = "h5ad/tahoe_log1p.h5ad"
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
    except Exception as e:
        debug("Failed to parse target_genes:", x, "err:", e)
        return []

drug_df["target_genes"] = drug_df["target_genes"].apply(parse_gene_list)

# -----------------------------
# Open h5ad
# -----------------------------
print("Opening h5ad in backed mode...")
adata = sc.read_h5ad(H5AD_PATH, backed="r")
debug("adata shape:", adata.shape)

var_tokens = adata.var.index.astype(int).values
debug("n_vars:", len(var_tokens))

# -----------------------------
# Core logic
# -----------------------------
results = []

for _, row in tqdm(drug_df.iterrows(), total=len(drug_df), desc="Processing drugs"):
    drug = row["drug"]
    target_genes = row["target_genes"]

    if not isinstance(target_genes, list) or len(target_genes) == 0:
        debug(drug, "skipped: empty target_genes")
        continue

    # map symbol -> token
    target_tokens = {}
    for g in target_genes:
        tok = symbol2token.get(g)
        if tok is None:
            debug(drug, "target gene not in gene_metadata:", g)
        else:
            target_tokens[g] = tok

    if not target_tokens:
        debug(drug, "skipped: no target genes mapped to token_id")
        continue

    # select all cells of this drug
    drug_mask = adata.obs["drug"] == drug
    n_drug_cells = int(drug_mask.sum())
    if n_drug_cells == 0:
        debug(drug, "skipped: no cells in h5ad")
        continue

    debug(drug, "cells:", n_drug_cells)
    drug_sub = adata[drug_mask].to_memory()

    for cell_line in drug_sub.obs["cell_line_id"].unique():
        cl_mask = drug_sub.obs["cell_line_id"] == cell_line
        n_cells = int(cl_mask.sum())
        if n_cells == 0:
            debug(drug, cell_line, "skipped: zero cells")
            continue

        sub = drug_sub[cl_mask]

        # delta per gene
        delta = sub.X.mean(axis=0) - sub.layers["ctrl_norm"].mean(axis=0)
        delta = np.asarray(delta).ravel()

        # rank by |delta|
        order = np.argsort(np.abs(delta))[::-1]
        ranked_tokens = var_tokens[order]
        token2rank = {tok: i + 1 for i, tok in enumerate(ranked_tokens)}

        # per-target-gene rank and top-K counts
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
            "n_cells": n_cells,
            "target_genes": ",".join(target_genes),
            "gene_ranks": gene_ranks,
        }
        # add top-K counts
        for k in TOPKS:
            record[f"n_in_top{k}"] = topk_counts[k]

        # best rank among all target genes
        record["best_rank"] = min(gene_ranks.values())

        results.append(record)

    del drug_sub

# -----------------------------
# Save results
# -----------------------------
res_df = pd.DataFrame(results)
res_df.to_csv("drug_target_in_topk_counts.csv", index=False)

print("\nDone! Total records:", len(res_df))
print(res_df.head())
