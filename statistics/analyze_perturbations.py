import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Analyze perturbation statistics")
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--gene-list", "-g", help="Optional: CSV or txt with 150 genes")
    parser.add_argument("--outdir", "-o", default="perturbation_analysis")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    print(f"读取数据：{args.input}")
    adata = sc.read_h5ad(args.input)
    print(adata)

    df = adata.obs[["target_gene", "batch"]].copy()
    df["is_control"] = df["target_gene"] == "non-targeting"

    if args.gene_list:
        target_genes = pd.read_csv(args.gene_list, header=None)[0].tolist()
        df = df[df["target_gene"].isin(target_genes)]

    print("\n=== 基因数量 ===")
    print(df["target_gene"].nunique())

    print("\n=== 每个 batch 的细胞数量 ===")
    batch_counts = df.groupby("batch").size().rename("cell_n").reset_index()

    print(batch_counts)

    # Count how many cells per batch per perturbation
    pert_batch_stats = df.groupby(["target_gene", "batch"]).size().rename("cell_n").reset_index()

    pert_stats = df.groupby("target_gene").size().rename("cell_n").reset_index()

    # ========== 计算扰动强度 ==========
    print("\n=== 计算扰动强度（表达变化基因数量） ===")

    # 基础表达矩阵
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    var_names = np.array(adata.var_names)

    control_mask = (adata.obs["target_gene"] == "non-targeting")
    X_ctrl = X[control_mask]

    results = []
    for gene in df["target_gene"].unique():
        if gene == "non-targeting": continue
        pert_mask = (adata.obs["target_gene"] == gene)
        X_pert = X[pert_mask]

        if X_pert.shape[0] < 5:  # avoid single-cell noise
            continue

        # 平均表达差值
        delta = np.mean(X_pert, axis=0) - np.mean(X_ctrl, axis=0)

        # threshold = 0.25 log-normalized expression
        significant = np.sum(np.abs(delta) > 0.25)

        if significant > 100:
            strength = "Strong"
        elif significant > 10:
            strength = "Medium"
        else:
            strength = "Weak"

        results.append({
            "gene": gene,
            "cell_count": X_pert.shape[0],
            "affected_genes": significant,
            "strength": strength
        })

    result_df = pd.DataFrame(results).sort_values("affected_genes", ascending=False)
    # result_df.to_csv(outdir / "perturbation_strength.csv", index=False)

    print("\n=== 扰动强度分类 ===")
    print(result_df["strength"].value_counts())

    # 保存
    batch_counts.to_csv(outdir / "batch_counts.csv", index=False)
    pert_batch_stats.to_csv(outdir / "perturbation_batch_counts.csv", index=False)
    pert_stats.to_csv(outdir / "gene_cell_counts.csv", index=False)

    print(f"\n文件已输出至 {outdir.resolve()}")

if __name__ == "__main__":
    main()
