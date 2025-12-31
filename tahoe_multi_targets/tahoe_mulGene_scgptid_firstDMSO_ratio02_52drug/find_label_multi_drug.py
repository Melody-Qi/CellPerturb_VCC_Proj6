import pandas as pd

INPUT_CSV = "all_conditions.csv"
OUT_CSV = "label_with_multiple_drugs.csv"

def main():
    df = pd.read_csv(INPUT_CSV)

    # groupby(label) + nunique(drug)
    label_drug_counts = (
        df.groupby("label")["drug"]
          .nunique()
          .reset_index(name="n_drugs")
    )

    # 只保留：同一个 label 对应多个 drug
    multi_drug_labels = label_drug_counts[
        label_drug_counts["n_drugs"] > 1
    ]

    # 把对应的具体 drug 一起列出来（方便人工检查）
    result = (
        df[df["label"].isin(multi_drug_labels["label"])]
        .groupby("label")
        .agg(
            drugs=("drug", lambda x: sorted(set(x))),
            n_drugs=("drug", "nunique"),
            genes=("target_gene", lambda x: sorted(set(x))),
            n_genes=("target_gene", "nunique"),
        )
        .reset_index()
    )

    # 保存
    result.to_csv(OUT_CSV, index=False)

    print("✅ Done")
    print(f"发现 {len(result)} 个 label 对应多个 drug")
    print(f"结果已保存到: {OUT_CSV}")

if __name__ == "__main__":
    main()
