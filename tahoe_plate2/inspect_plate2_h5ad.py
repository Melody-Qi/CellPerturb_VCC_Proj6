import os
import scanpy as sc
import pandas as pd

H5AD_PATH = "plate2_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
SAMPLE_N = 1000   # 安全抽样，不会炸内存


def human_readable_size(path):
    size = os.path.getsize(path)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} PB"


def main():
    print("=" * 80)
    print("Inspecting h5ad file:")
    print(H5AD_PATH)

    # --------------------------------------------------
    # 1. 文件大小
    # --------------------------------------------------
    if not os.path.exists(H5AD_PATH):
        raise FileNotFoundError(H5AD_PATH)

    print("\n[1] File size")
    print("    ", human_readable_size(H5AD_PATH))

    # --------------------------------------------------
    # 2. 只加载结构（backed 模式）
    # --------------------------------------------------
    print("\n[2] Loading AnnData structure (backed='r') ...")
    adata = sc.read_h5ad(H5AD_PATH, backed="r")

    print("\n[3] Basic AnnData info")
    print(f"    n_cells (obs): {adata.n_obs}")
    print(f"    n_genes (var): {adata.n_vars}")

    # --------------------------------------------------
    # 3. obs / var / uns / layers
    # --------------------------------------------------
    print("\n[4] obs columns")
    print(list(adata.obs.keys()))

    print("\n[5] var columns")
    print(list(adata.var.keys()))

    print("\n[6] layers")
    print(list(adata.layers.keys()))

    print("\n[7] uns keys")
    print(list(adata.uns.keys()))

    # --------------------------------------------------
    # 4. 关键字段分布（不触碰 X）
    # --------------------------------------------------
    def safe_value_counts(series, name, topk=10):
        print(f"\n[Distribution] {name}")
        if name not in adata.obs:
            print("  -> NOT FOUND")
            return
        vc = adata.obs[name].value_counts().head(topk)
        print(vc)

    safe_value_counts("drug")
    safe_value_counts("cell_line")
    safe_value_counts("plate")

    # --------------------------------------------------
    # 5. 安全抽样一小块到内存
    # --------------------------------------------------
    print("\n[8] Loading small subset into memory")
    n = min(SAMPLE_N, adata.n_obs)
    adata_small = adata[:n].to_memory()

    print("    subset shape:", adata_small.shape)
    print("    X type:", type(adata_small.X))

    # 查看前几条样本
    print("\n[9] First 5 obs rows (subset)")
    print(adata_small.obs.head())

    # --------------------------------------------------
    # 6. 是否包含 ADAGRASIB（示例）
    # --------------------------------------------------
    if "drug" in adata_small.obs:
        n_ada = (adata_small.obs["drug"] == "ADAGRASIB").sum()
        print(f"\n[10] ADAGRASIB in first {n} cells: {n_ada}")

    print("\nDone. Structure inspection finished.")
    print("=" * 80)


if __name__ == "__main__":
    main()
