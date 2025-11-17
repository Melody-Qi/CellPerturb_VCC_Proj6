#!/usr/bin/env python3
"""
compute_perturbation_strength.py

Produces perturbation_strength.csv with columns:
gene,cell_count,affected_genes,strength

Usage examples:
# default: auto top-150 genes from adata (exclude non-targeting)
python compute_perturbation_strength.py -i ../STATE/preprocessed_data/preprocessed_training_data_2000.h5ad -o ./pert_strength_out

# provide your own gene list (one gene per line)
python compute_perturbation_strength.py -i ../STATE/... -g my150genes.txt -o ./pert_strength_out
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")


def bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment"""
    _, p_adj, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    return p_adj


def ensure_matrix(adata):
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def compute_strength(adata, gene_list, min_cells=20, logfc_thresh=0.25):
    X = ensure_matrix(adata)
    obs = adata.obs
    genes_all = np.array(adata.var_names)

    # control mask
    if 'target_gene' not in obs.columns:
        raise RuntimeError("adata.obs must contain 'target_gene' column (control labeled 'non-targeting')")

    ctrl_mask = (obs['target_gene'] == 'non-targeting').values
    n_ctrl = int(ctrl_mask.sum())
    if n_ctrl < 5:
        raise RuntimeError("Too few control cells for DE (found {}).".format(n_ctrl))

    ctrl_mean = X[ctrl_mask].mean(axis=0)

    results = []
    for g in gene_list:
        # boolean mask of cells for this perturbation
        pert_mask = (obs['target_gene'] == g).values
        n_cells = int(np.sum(pert_mask))
        if n_cells < min_cells:
            results.append((g, n_cells, 0, 'TooFewCells'))
            continue

        # means and logFC (on whatever scale X is: if already log1p, then logFC is difference)
        pert_mean = X[pert_mask].mean(axis=0)
        logfc = pert_mean - ctrl_mean  # per-gene

        # t-test across genes (vectorized)
        tstat, pvals = stats.ttest_ind(X[pert_mask], X[ctrl_mask], axis=0, equal_var=False, nan_policy='omit')
        pvals = np.nan_to_num(pvals, nan=1.0, posinf=1.0, neginf=1.0)

        # FDR adjust
        p_adj = bh_adjust(pvals)

        # significant genes
        sig_mask = (p_adj < 0.05) & (np.abs(logfc) >= logfc_thresh)
        n_sig = int(np.sum(sig_mask))

        # append: gene, cell_count, affected_genes, strength
        if n_sig > 100:
            strength = 'Strong'
        elif n_sig >= 10:
            strength = 'Medium'
        else:
            strength = 'Weak'

        results.append((g, n_cells, n_sig, strength))

    df = pd.DataFrame(results, columns=['gene','cell_count','affected_genes','strength'])
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", required=True, help="path to h5ad")
    parser.add_argument("--genes", "-g", default=None, help="optional gene list (one per line). If not provided, top 150 pert genes by frequency are used.")
    parser.add_argument("--outdir", "-o", default=".", help="output directory")
    parser.add_argument("--min_cells", type=int, default=20, help="min pert cells to run DE")
    parser.add_argument("--logfc", type=float, default=0.25, help="logFC threshold for calling DE (same scale as X)")
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("读取数据：", str(input_path))
    adata = sc.read_h5ad(str(input_path))
    print(adata)

    # decide gene list
    if args.genes:
        gene_list = [g.strip() for g in Path(args.genes).read_text().splitlines() if g.strip()]
    else:
        counts = adata.obs['target_gene'].value_counts()
        if 'non-targeting' in counts.index:
            counts = counts.drop('non-targeting')
        topn = 150
        gene_list = list(counts.index[:topn].astype(str))
        print(f"自动选取 top-{topn} perturbations (by frequency).")

    print(f"计算 {len(gene_list)} perturbations ... (min_cells={args.min_cells}, logfc_thresh={args.logfc})")
    df = compute_strength(adata, gene_list, min_cells=args.min_cells, logfc_thresh=args.logfc)

    out_csv = outdir / "perturbation_strength.csv"
    df.to_csv(out_csv, index=False)
    print("已写出：", out_csv.resolve())


if __name__ == "__main__":
    main()
