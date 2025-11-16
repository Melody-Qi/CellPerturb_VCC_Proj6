import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Training data statistics")
    parser.add_argument("--input", "-i", required=True, help="Path to the h5ad file")
    parser.add_argument("--outdir", "-o", default=".", help="Output directory")
    args = parser.parse_args()

    input_path = Path(args.input)
    prefix = input_path.stem  # e.g. "preprocessed_training_data_2000"
    
    # 输出目录结构: <dataset_name>/
    output_dir = Path(args.outdir) / prefix
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    print(f"读取数据：{args.input}")

    adata = sc.read_h5ad(args.input)

    print(adata)
    print("obs keys:", adata.obs.keys())
    print("var keys:", adata.var.keys())

    # X
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

    # Stats
    gene_means = np.mean(X, axis=0)
    gene_max = np.max(X, axis=0)
    gene_min = np.min(X, axis=0)
    gene_nonzero_fraction = np.sum(X > 0, axis=0) / X.shape[0]

    print("平均表达值均值：", np.mean(gene_means))
    print("平均最大表达值：", np.mean(gene_max))
    print("平均最小值表达值：", np.mean(gene_min))
    print("基因完全不表达的比例：", np.mean(gene_nonzero_fraction == 0))

    # Save CSV with prefix
    out_csv = output_dir / f"{prefix}_gene_statistics.csv"
    df_stats = pd.DataFrame({
        "gene": adata.var_names,
        "mean": gene_means,
        "max": gene_max,
        "min": gene_min,
        "nonzero_fraction": gene_nonzero_fraction
    })
    df_stats.to_csv(out_csv, index=False)
    print(f"保存统计文件: {out_csv}")

    # Plot
    plt.figure(figsize=(8,5))
    plt.hist(np.log1p(gene_means), bins=50)
    plt.xlabel("log(1 + gene_means)")
    plt.ylabel("gene_num")
    plt.title(f"{len(adata.var_names)} genes mean distribution")
    out_plot = plot_dir / f"{prefix}_gene_mean_distribution.png"
    plt.savefig(out_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"保存统计图: {out_plot}")

if __name__ == "__main__":
    main()
