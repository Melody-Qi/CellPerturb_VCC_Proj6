import pandas as pd

# 读取 sample metadata
df = pd.read_parquet("sample_metadata.parquet")

print("Total samples:", len(df))
print("Unique drugs:", df["drug"].nunique())

# 基础清洗
df["drug_clean"] = (
    df["drug"]
    .astype(str)
    .str.strip()
    .str.upper()
)

# 去重
unique_drugs = (
    df[["drug_clean"]]
    .drop_duplicates()
    .sort_values("drug_clean")
    .reset_index(drop=True)
)

print("Unique cleaned drugs:", len(unique_drugs))
unique_drugs.head()

unique_drugs.to_csv(
    "tahoe_drugs.csv",
    index=False
)

print("Saved tahoe_drugs.csv")

