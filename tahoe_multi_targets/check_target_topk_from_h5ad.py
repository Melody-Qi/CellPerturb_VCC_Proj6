#!/usr/bin/env python
"""
Evaluate whether drug target genes (single or multi-target)
appear in top-K differential expression ranks from Tahoe h5ad.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from tqdm import tqdm


H5AD_PATH = "tahoe_mulGene_scgptid_firstDMSO_ratio02_52drug/h5ad/tahoe_log1p.h5ad"
OUT_CSV = "target_topk_from_h5ad.csv"
TOPKS = [5, 10, 50, 100]


def split_targets(target_str: str) -> list[str]:
    """Split canonical target_gene string into list."""
    if target_str in ("", "ctrl", None):
        return []
    return [g.strip() for g in str(target_str).split("+") if g.strip()]


def main():
    print(f"Loading {H5AD_PATH}")
    adata = sc.read_h5ad(H5AD_PATH)

    # gene token ids (string!)
    gene_ids = adata.var.index.to_numpy()

    # build fast index lookup
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}

    results = []

    # group by condition (cell_line × drug)
    grouped = adata.obs.groupby(["cell_line_id", "drug", "target_gene"])

    for (cell, drug, target_gene), obs_idx in tqdm(
        grouped.indices.items(),
        desc="Processing conditions"
    ):
        if drug == "DMSO_TF":
            continue

        targets = split_targets(target_gene)
        if not targets:
            continue

        # map target gene → var index
        target_indices = []
        for t in targets:
            if t in gene_to_idx:
                target_indices.append(gene_to_idx[t])

        if not target_indices:
            continue

        sub = adata[obs_idx]

        # mean perturbation
        X = sub.X
        if issparse(X):
            pert_mean = np.asarray(X.mean(axis=0)).ravel()
        else:
            pert_mean = X.mean(axis=0)

        # mean control
        ctrl = sub.layers["ctrl_norm"]
        if issparse(ctrl):
            ctrl_mean = np.asarray(ctrl.mean(axis=0)).ravel()
        else:
            ctrl_mean = ctrl.mean(axis=0)

        delta = pert_mean - ctrl_mean
        ranks = np.argsort(-np.abs(delta))  # descending

        # compute rank for each target
        target_ranks = {
            gene_ids[i]: int(np.where(ranks == i)[0][0]) + 1
            for i in target_indices
        }

        row = {
            "cell_line_id": cell,
            "drug": drug,
            "target_gene": target_gene,
            "n_targets": len(target_indices),
            "n_cells": len(obs_idx),
        }

        for k in TOPKS:
            row[f"top{k}"] = any(r <= k for r in target_ranks.values())

        # optional: store min rank
        row["best_rank"] = min(target_ranks.values())

        results.append(row)

    df = pd.DataFrame(results)
    df.to_csv(OUT_CSV, index=False)

    print(f"\n✅ Done. Results saved to {OUT_CSV}")
    print("\nSummary:")
    for k in TOPKS:
        print(f"Top{k} hit rate: {df[f'top{k}'].mean():.3f}")


if __name__ == "__main__":
    main()
