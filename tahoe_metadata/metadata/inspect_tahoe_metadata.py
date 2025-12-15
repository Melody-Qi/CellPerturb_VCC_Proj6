import pandas as pd

# =========================
# 1. 读取 metadata
# =========================

sample_md = pd.read_parquet(
    "sample_metadata.parquet"
)

obs_md = pd.read_parquet(
    "obs_metadata.parquet"
)

print("=== Loaded metadata ===")
print(f"Sample metadata shape: {sample_md.shape}")
print(f"Obs (cell) metadata shape: {obs_md.shape}")

# =========================
# 2. 查看 sample_metadata
# =========================

print("\n=== Sample metadata columns ===")
print(sample_md.columns.tolist())

print("\n=== Sample metadata preview ===")
print(sample_md.head())

print("\nNumber of samples:", len(sample_md))
print("Number of unique drugs:", sample_md["drug"].nunique())

print("\nTop drugs by sample count:")
print(sample_md["drug"].value_counts().head(20))

# =========================
# 3. 查看 obs_metadata
# =========================

print("\n=== Obs metadata columns ===")
print(obs_md.columns.tolist())

print("\n=== Obs metadata preview ===")
print(obs_md.head())

print("\nNumber of cells:", len(obs_md))
print("Cell lines:")
print(obs_md["cell_line"].value_counts())

# =========================
# 4. sample ↔ cell 的对应关系
# =========================

print("\n=== How many cells per sample? ===")
cells_per_sample = obs_md["sample"].value_counts()
print(cells_per_sample.head())

# =========================
# 5. 示例：选一个 drug，看它有哪些 sample / cell_line
# =========================

example_drug = sample_md["drug"].value_counts().index[1]
print(f"\n=== Example drug: {example_drug} ===")

drug_samples = sample_md[sample_md["drug"] == example_drug]
print("Samples for this drug:")
print(drug_samples[["sample", "drug", "plate"]])

drug_cells = obs_md[obs_md["sample"].isin(drug_samples["sample"])]
print("\nCell lines for this drug:")
print(drug_cells["cell_line"].value_counts())

print("\nDone.")
