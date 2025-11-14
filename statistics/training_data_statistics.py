import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description="Training data statistics")
    parser.add_argument("--input", "-i", required=True, help="Path to the h5ad file")
    args = parser.parse_args()

    print(f"读取数据：{args.input}")

    adata = sc.read_h5ad(args.input)

    # 查看基本信息
    print(adata)
    print("obs keys:", adata.obs.keys())
    print("var keys:", adata.var.keys())

    # 表达矩阵
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

    # 基本统计
    gene_means = np.mean(X, axis=0)
    gene_max = np.max(X, axis=0)
    gene_min = np.min(X, axis=0)
    gene_nonzero_fraction = np.sum(X > 0, axis=0) / X.shape[0]

    print("平均表达值均值：", np.mean(gene_means))
    print("平均最大表达值：", np.mean(gene_max))
    print("平均最小值表达值：", np.mean(gene_min))
    print("基因完全不表达的比例：", np.mean(gene_nonzero_fraction == 0))

    df_stats = pd.DataFrame({
    "gene": adata.var_names,
    "mean": gene_means,
    "max": gene_max,
    "min": gene_min,
    "nonzero_fraction": gene_nonzero_fraction
    })
    df_stats.to_csv("gene_statistics.csv", index=False)


    import os

    # 创建保存目录
    os.makedirs("plots", exist_ok=True)

    plt.figure(figsize=(8,5))
    plt.hist(np.log1p(gene_means), bins=50)
    plt.xlabel("log(1 + gene_means)")
    plt.ylabel("gene_num")
    plt.title("18080 genes mean distribution")

    # 保存图像
    plt.savefig("plots/gene_mean_distribution.png", dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()
