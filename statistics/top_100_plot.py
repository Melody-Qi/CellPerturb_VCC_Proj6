import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# User configuration
# ----------------------------

HVG_ONLY = True
TOP_K = 100
SAVE_DIR = "./perturbation_analysis/"

# ----------------------------
# Load data
# ----------------------------
adata = sc.read_h5ad("../STATE/preprocessed_data/preprocessed_training_data_2000.h5ad")
print("Loaded:", adata)

# Identify control cells
ctrl_mask = (adata.obs["target_gene"] == "non-targeting")

pert_genes = [g for g in adata.obs["target_gene"].unique() if g != "non-targeting"]

# HVG mask
if HVG_ONLY and "highly_variable" in adata.var:
    hvg_mask = adata.var["highly_variable"].values
else:
    hvg_mask = np.ones(adata.n_vars, dtype=bool)

os.makedirs(SAVE_DIR, exist_ok=True)

# global coordinate limit
global_max = float(adata.X.max())
print("global max =", global_max)

# precompute control mean
ctrl_mean = np.array(adata[ctrl_mask].X.mean(axis=0)).ravel()

# ----------------------------
# Process each perturbation
# ----------------------------
for pert in pert_genes:
    print(f"\n=== Processing {pert} ===")

    sub_mask = (adata.obs["target_gene"] == pert)
    sub = adata[sub_mask]

    pert_mean = np.array(sub.X.mean(axis=0)).ravel()
    delta = pert_mean - ctrl_mean

    # top-100
    top_gene_idx = np.argsort(np.abs(delta))[::-1][:TOP_K]

    # output folder
    out_dir = os.path.join(SAVE_DIR, pert)
    os.makedirs(out_dir, exist_ok=True)

    # --------------------------------------------------
    # 1) Save full 18,080-gene statistics (sorted by |delta|)
    # --------------------------------------------------
    full_df = pd.DataFrame({
        "gene": adata.var_names,
        "ctrl_mean": ctrl_mean,
        "pert_mean": pert_mean,
        "delta": delta,
        "is_HVG": hvg_mask,
    })

    # compute variance
    full_df["ctrl_var"] = adata[ctrl_mask].X.toarray().var(axis=0)
    full_df["pert_var"] = sub.X.toarray().var(axis=0)

    # sort by |delta|: descending
    full_df = full_df.iloc[np.argsort(np.abs(delta))[::-1]]

    full_df.to_csv(os.path.join(out_dir, f"{pert}_18080_stats.csv"), index=False)

    # --------------------------------------------------
    # 2) Plot top-100 genes
    # --------------------------------------------------
    for gi in top_gene_idx:
        gene = adata.var_names[gi]
        tag = "HVG" if hvg_mask[gi] else ""

        ctrl_vals = adata[ctrl_mask].X[:, gi].toarray().ravel()
        pert_vals = sub.X[:, gi].toarray().ravel()

        plt.figure(figsize=(5, 4))
        # plt.scatter(ctrl_vals, pert_vals, s=5, alpha=0.5)

        plt.violinplot([ctrl_vals, pert_vals], showmeans=True)

        plt.xticks([1, 2], ["Control", "Perturbed"])
        plt.ylabel("Expression")
        plt.title(f"{pert} — {gene} {tag}")

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{gene}_violin.png"), dpi=120)
        plt.close()

print("\nAll tasks complete.")
