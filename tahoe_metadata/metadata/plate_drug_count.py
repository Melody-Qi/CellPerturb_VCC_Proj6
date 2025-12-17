import pandas as pd

df = pd.read_parquet("sample_metadata.parquet")

plate_drug_count = (
    df.groupby("plate")["drug"]
    .nunique()
    .sort_values(ascending=False)
)

print(plate_drug_count)
