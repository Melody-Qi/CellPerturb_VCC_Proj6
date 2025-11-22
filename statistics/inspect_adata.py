import anndata as ad
import pandas as pd

def inspect_adata(path):
    print(f"📥 Loading AnnData from: {path}")
    adata = ad.read_h5ad(path)

    print("\n=== 📌 Basic Info ===")
    print(adata)

    print("\n=== 🔑 obs columns ===")
    print(list(adata.obs.columns))

    print("\n=== 🔍 obs (前 10 行) ===")
    print(adata.obs.head(10))

    print("\n=== 🔑 var columns ===")
    print(list(adata.var.columns))

    print("\n=== 🔍 var (前 10 行) ===")
    print(adata.var.head(10))

    print("\n=== 📦 X matrix type & shape ===")
    print(type(adata.X), adata.X.shape)

    print("\n=== 🧬 uns keys ===")
    print(list(adata.uns.keys()))

    print("\n=== FULL OBS TABLE INFO ===")
    print(adata.obs.info())

    print("\n=== OBS HEAD (no truncated columns) ===")
    print(pd.DataFrame(adata.obs).head(10).to_string())

    print("\n=== VAR INFO ===")
    print(adata.var.info())

    print("\n=== VAR HEAD (should show only index) ===")
    print(pd.DataFrame(adata.var).head(10).to_string())

if __name__ == "__main__":
    inspect_adata("../STATE/competition_support_set/competition_train.h5")
