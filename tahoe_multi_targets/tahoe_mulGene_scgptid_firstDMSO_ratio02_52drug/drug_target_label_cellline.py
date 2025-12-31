import glob
import json
import pandas as pd

PARQUET_FILES = sorted(glob.glob("single_target_*.parquet"))
OUT_CSV = "all_conditions.csv"

# ===============================
# global set for unique conditions
# ===============================
unique_conditions = set()

# ===============================
# main loop
# ===============================
for f in PARQUET_FILES:
    print(f"Processing {f}")

    df = pd.read_parquet(
        f,
        columns=["drug", "target_gene", "label", "cell_line_id"]
    )

    for _, row in df.iterrows():
        drug = row["drug"]
        target_gene = row["target_gene"]
        label = int(row["label"])
        cell = row["cell_line_id"]

        unique_conditions.add((drug, target_gene, label, cell))

    del df

# ===============================
# write to CSV
# ===============================
out_df = pd.DataFrame(
    list(unique_conditions),
    columns=["drug", "target_gene", "label", "cell_line_id"]
)

out_df.sort_values(
    ["drug", "cell_line_id", "target_gene"],
    inplace=True
)

out_df.to_csv(OUT_CSV, index=False)

print(f"\n✅ Done! {len(out_df)} unique conditions written to {OUT_CSV}")
