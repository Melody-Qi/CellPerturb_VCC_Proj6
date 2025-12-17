import pandas as pd

# =========================
# 1. 读取 metadata
# =========================
sample_md = pd.read_parquet("sample_metadata.parquet")
obs_md = pd.read_parquet("obs_metadata.parquet")

PLATE = "plate2"

print(f"=== Inspecting {PLATE} ===")

# =========================
# 2. plate2 的 sample 级信息
# =========================
plate2_samples = sample_md[sample_md["plate"] == PLATE]

print("\n[1] Samples in plate2")
print(f"Number of samples: {len(plate2_samples)}")

print("\nSample → drug / drugname_drugconc")
print(
    plate2_samples[
        ["sample", "drug", "drugname_drugconc"]
    ].sort_values("drug")
)

# =========================
# 3. plate2 的 cell 级信息
# =========================
plate2_cells = obs_md[obs_md["plate"] == PLATE]

print("\n[2] Cells in plate2")
print(f"Number of cells: {len(plate2_cells)}")

# =========================
# 4. drug → sample → cell_line 汇总表
# =========================
summary = (
    plate2_cells
    .groupby(["drug", "drugname_drugconc", "sample", "cell_line"])
    .size()
    .reset_index(name="n_cells")
    .sort_values(["drug", "cell_line"], ascending=[True, False])
)

print("\n[3] drug → sample → cell_line (with cell counts)")
print(summary)

# =========================
# 5. 每个 drug 覆盖哪些 cell_line
# =========================
drug_cellline = (
    plate2_cells
    .groupby("drug")["cell_line"]
    .unique()
    .reset_index()
)

print("\n[4] drug → cell_lines")
print(drug_cellline)

# =========================
# 6. 保存结果（强烈建议）
# =========================
summary.to_csv("plate2_drug_sample_cellline_summary.csv", index=False)
drug_cellline.to_csv("plate2_drug_celllines.csv", index=False)

print("\nSaved:")
print(" - plate2_drug_sample_cellline_summary.csv")
print(" - plate2_drug_celllines.csv")
print("\nDone.")
