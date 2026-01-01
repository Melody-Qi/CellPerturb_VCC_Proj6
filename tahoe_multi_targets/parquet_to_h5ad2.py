#!/usr/bin/env python
"""
Safe conversion from Tahoe single_target_*.parquet to per-file H5AD.

Key features:
- CSR construction (no LIL)
- Batchable via --start / --end
- OS-safe (no OOM killer)
"""

import argparse
import glob
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from tqdm import tqdm


# -----------------------------
# Args
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    return parser.parse_args()


# -----------------------------
# Build global gene vocabulary
# -----------------------------
def build_gene_vocabulary(files):
    all_genes = set()
    for f in tqdm(files, desc="Building gene vocabulary"):
        df = pd.read_parquet(f, columns=["genes"])
        for g in df["genes"]:
            all_genes.update(g)
        del df
    sorted_genes = sorted(all_genes)
    return {g: i for i, g in enumerate(sorted_genes)}


# -----------------------------
# Process one parquet → H5AD
# -----------------------------
def process_one_parquet(filepath, gene_to_idx, n_genes, out_dir):
    df = pd.read_parquet(filepath)
    n_cells = len(df)

    rows_x, cols_x, data_x = [], [], []
    rows_c, cols_c, data_c = [], [], []

    for i in range(n_cells):
        row = df.iloc[i]

        # treated
        for g, e in zip(row["genes"], row["expressions"]):
            idx = gene_to_idx.get(g)
            if idx is not None:
                rows_x.append(i)
                cols_x.append(idx)
                data_x.append(e)

        # control
        for g, e in zip(row["ctrl_genes"], row["ctrl_expressions"]):
            idx = gene_to_idx.get(g)
            if idx is not None:
                rows_c.append(i)
                cols_c.append(idx)
                data_c.append(e)

    X = csr_matrix((data_x, (rows_x, cols_x)), shape=(n_cells, n_genes))
    ctrl = csr_matrix((data_c, (rows_c, cols_c)), shape=(n_cells, n_genes))

    obs = pd.DataFrame(
        {
            "condition": df["target_gene"] + "+" + df["drug"],
            "cell_line_id": df["cell_line_id"].astype(str),
            "plate": df["plate"].astype(str),
            "sample": df["sample"].astype(str),
            "label": df["label"],
            "overlap_ratio": df["overlap_ratio"],
            "control": 0,
            "target_gene": df["target_gene"].astype(str),
            "drug": df["drug"].astype(str),
            "ctrl_drug": df["ctrl_drug"].astype(str),
        }
    )
    obs["condition_name"] = obs["cell_line_id"] + "_" + obs["condition"]

    var = pd.DataFrame(index=[str(i) for i in range(n_genes)])
    var["gene_idx"] = list(range(n_genes))

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=var,
        layers={"ctrl": ctrl},
    )

    for col in obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")

    out_path = out_dir / f"{Path(filepath).stem}.h5ad"
    adata.write_h5ad(out_path)

    # hard cleanup
    del df, X, ctrl, adata, rows_x, cols_x, data_x, rows_c, cols_c, data_c


# -----------------------------
# Main
# -----------------------------
def main():
    args = parse_args()

    data_dir = Path("tahoe_mulGene_scgptid_firstDMSO_ratio02_52drug")
    out_dir = data_dir / "h5ad"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_files = sorted(glob.glob(str(data_dir / "single_target_*.parquet")))
    files = all_files[args.start : args.end]

    print(f"Processing files {args.start}:{args.end} ({len(files)} files)")

    print("Building gene vocabulary (global)...")
    gene_to_idx = build_gene_vocabulary(all_files)
    n_genes = len(gene_to_idx)
    print(f"Gene vocab size: {n_genes}")

    for f in tqdm(files, desc="Converting to h5ad"):
        process_one_parquet(f, gene_to_idx, n_genes, out_dir)

    print("Batch done.")


if __name__ == "__main__":
    main()
