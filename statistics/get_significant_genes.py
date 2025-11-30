import pandas as pd
import os

# 1. 读取 DE 结果表
df = pd.read_csv("../vcc_data/cell-eval-outdir/real_de.csv")

# 2. 筛选显著基因（FDR ≤ 0.05）
df_sig = df[df["fdr"] <= 0.05]

# 3. 总输出目录
base_dir = "./perturbation_analysis"
os.makedirs(base_dir, exist_ok=True)

# 4. 按 target 分组并保存
for target, sub_df in df_sig.groupby("target"):
    # 为每个 target 创建子文件夹
    target_dir = os.path.join(base_dir, target)
    os.makedirs(target_dir, exist_ok=True)
    
    # 保存文件
    out_path = os.path.join(target_dir, f"{target}.csv")
    sub_df.to_csv(out_path, index=False)
    
    print(f"✔ Saved {out_path} (n={len(sub_df)})")

print("All done.")
