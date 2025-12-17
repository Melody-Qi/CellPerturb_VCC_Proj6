import pandas as pd

# 读两个 metadata
sample_df = pd.read_parquet("sample_metadata.parquet")
obs_df = pd.read_parquet("obs_metadata.parquet")

# 清洗 drug
sample_df["drug_clean"] = sample_df["drug"].str.upper().str.strip()

target_drug = "ADAGRASIB"

# 1️⃣ 找 ADAGRASIB 对应的 sample
target_samples = sample_df.loc[
    sample_df["drug_clean"] == target_drug, "sample"
].unique()

print(target_samples)

print(f"Number of samples with {target_drug}: {len(target_samples)}")

# 2️⃣ 在 obs_metadata 里筛这些 sample
obs_hits = obs_df[obs_df["sample"].isin(target_samples)]

# 3️⃣ 看 cell line 分布
cell_line_counts = (
    obs_hits.groupby(["plate", "cell_line"])
    .size()
    .sort_values(ascending=False)
)

print(cell_line_counts.head(20))
