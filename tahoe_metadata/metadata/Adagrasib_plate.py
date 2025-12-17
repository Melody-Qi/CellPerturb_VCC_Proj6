import pandas as pd

# 读取 sample-level metadata
df = pd.read_parquet("sample_metadata.parquet")

# 清洗 drug 名称
df["drug_clean"] = (
    df["drug"]
    .astype(str)
    .str.strip()
    .str.upper()
)

target_drug = "ADAGRASIB"

hits = df[df["drug_clean"] == target_drug]

print(f"Total samples with {target_drug}: {len(hits)}")

plates = (
    hits.groupby("plate")["sample"]
    .nunique()
    .sort_values(ascending=False)
)

print("\nPlates containing ADAGRASIB:")
print(plates)
