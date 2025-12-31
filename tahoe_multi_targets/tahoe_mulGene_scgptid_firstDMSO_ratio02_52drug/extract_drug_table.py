import pandas as pd

INPUT_CSV = "all_conditions.csv"
OUT_CSV = "drug_summary.csv"

def uniq_sorted(x):
    """去重 + 排序，保证稳定可复现"""
    return sorted(set(x))

def main():
    df = pd.read_csv(INPUT_CSV)

    # 必要列检查（防止 silent bug）
    required_cols = {"drug", "target_gene", "label", "cell_line_id"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # 以 drug 为 key 聚合
    drug_df = (
        df.groupby("drug")
          .agg(
              target_genes=("target_gene", uniq_sorted),
              labels=("label", uniq_sorted),
              cell_line_ids=("cell_line_id", uniq_sorted),
              n_targets=("target_gene", "nunique"),
              n_labels=("label", "nunique"),
              n_cell_lines=("cell_line_id", "nunique"),
          )
          .reset_index()
    )

    drug_df.to_csv(OUT_CSV, index=False)

    print("✅ Done")
    print(f"共提取 {len(drug_df)} 个不重复 drug")
    print(f"结果已保存到: {OUT_CSV}")

if __name__ == "__main__":
    main()
