#!/usr/bin/env python3
"""
statistics/gene_delta_statistics.py

用法（在 statistics/ 目录下运行）：
    python gene_delta_statistics.py -i ../STATE/vcc_data/adata_Training.h5ad --topk 50

输出（写到当前目录 statistics/):
  - gene_delta_all_sorted.csv
  - gene_delta_top50.csv
  - gene_delta_with_nonzero_fraction.csv
  - plots/
      - gene_delta_distribution.png
      - ctrl_vs_pert_scatter.png
      - top50_delta_barh.png
"""
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def safe_toarray(X):
    try:
        if hasattr(X, "toarray"):
            return X.toarray()
        else:
            return np.asarray(X)
    except MemoryError:
        raise MemoryError("数据太大，无法转换为 dense array。请在有足够内存的机器上运行，或改用稀疏运算。")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", "-i", required=True, help="输入 h5ad 文件路径（相对或绝对）")
    p.add_argument("--topk", type=int, default=50, help="保存 top K 的基因数量（默认50）")
    args = p.parse_args()

    base_dir = os.getcwd()   # 默认输出到当前工作目录（你在 statistics/ 下运行）
    plots_dir = os.path.join(base_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("读取数据：", args.input)
    adata = sc.read_h5ad(args.input)
    print(adata)
    if "target_gene" not in adata.obs.columns:
        raise KeyError("adata.obs 中没有 'target_gene' 列，请确认 h5ad 的列名。")

    ctrl = adata[adata.obs["target_gene"] == "non-targeting"]
    pert = adata[adata.obs["target_gene"] != "non-targeting"]
    print(f"总细胞: {adata.n_obs}, control: {ctrl.n_obs}, perturbed: {pert.n_obs}")

    X_ctrl = safe_toarray(ctrl.X)
    X_pert = safe_toarray(pert.X)

    ctrl_means = np.mean(X_ctrl, axis=0)
    pert_means = np.mean(X_pert, axis=0)
    delta = pert_means - ctrl_means
    abs_delta = np.abs(delta)

    df = pd.DataFrame({
        "gene": np.array(adata.var_names, dtype=str),
        "ctrl_mean": ctrl_means,
        "pert_mean": pert_means,
        "delta": delta,
        "abs_delta": abs_delta
    })
    df_sorted = df.sort_values("abs_delta", ascending=False).reset_index(drop=True)

    all_csv = os.path.join(base_dir, "gene_delta_all_sorted.csv")
    df_sorted.to_csv(all_csv, index=False)
    print("已保存：", all_csv)

    topk = args.topk
    df_topk = df_sorted.head(topk)
    topk_csv = os.path.join(base_dir, f"gene_delta_top{topk}.csv")
    df_topk.to_csv(topk_csv, index=False)
    print(f"已保存 top{topk}：", topk_csv)
    print("Top genes:", list(df_topk["gene"].values))

    ctrl_nonzero_fraction = np.sum(X_ctrl > 0, axis=0) / max(1, X_ctrl.shape[0])

    never_expressed_genes = adata.var_names[ctrl_nonzero_fraction == 0]
    print("Control 中永不表达的基因数量：", len(never_expressed_genes))

    df_extra = df_sorted.copy()
    df_extra["ctrl_nonzero_fraction"] = ctrl_nonzero_fraction[df_extra.index]
    df_nonzero_csv = os.path.join(base_dir, "gene_delta_with_nonzero_fraction.csv")
    df_extra.to_csv(df_nonzero_csv, index=False)
    print("已保存（包含 ctrl 非零比例）：", df_nonzero_csv)

    # 绘图 1: delta 分布
    plt.figure(figsize=(8,5))
    plt.hist(df_sorted["delta"], bins=120)
    plt.xlabel("pert_mean - ctrl_mean")
    plt.ylabel("gene count")
    plt.title("Gene expression change distribution (pert - ctrl)")
    f1 = os.path.join(plots_dir, "gene_delta_distribution.png")
    plt.savefig(f1, dpi=300, bbox_inches="tight")
    plt.close()
    print("已保存：", f1)

    # 绘图 2: 控制 vs 扰动 平均表达散点（log1p）
    plt.figure(figsize=(6,6))
    plt.scatter(np.log1p(df_sorted["ctrl_mean"]), np.log1p(df_sorted["pert_mean"]), s=6, alpha=0.5)
    plt.xlabel("log1p(ctrl_mean)")
    plt.ylabel("log1p(pert_mean)")
    plt.title("Control vs Perturbed Mean (per gene)")
    f2 = os.path.join(plots_dir, "ctrl_vs_pert_scatter.png")
    plt.savefig(f2, dpi=300, bbox_inches="tight")
    plt.close()
    print("已保存：", f2)

    # 绘图 3: topK 条形图（delta）
    plt.figure(figsize=(10,6))
    plt.barh(df_topk["gene"][::-1], df_topk["delta"][::-1])
    plt.xlabel("pert_mean - ctrl_mean")
    plt.title(f"Top {topk} genes by abs_delta")
    f3 = os.path.join(plots_dir, f"top{topk}_delta_barh.png")
    plt.savefig(f3, dpi=300, bbox_inches="tight")
    plt.close()
    print("已保存：", f3)

    print("全部完成，输出在当前目录（base）:", base_dir)
    print("图片在：", plots_dir)

if __name__ == "__main__":
    main()
