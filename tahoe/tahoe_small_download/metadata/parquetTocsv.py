import pandas as pd
import os

# 假设你的 parquet 文件路径
input_dir = "."  # parquet 所在文件夹
output_dir = "./csv_output"  # 输出 csv 的文件夹

os.makedirs(output_dir, exist_ok=True)

# 要转换的文件列表
parquet_files = ["drug_metadata.parquet", "gene_metadata.parquet"]

for pq_file in parquet_files:
    pq_path = os.path.join(input_dir, pq_file)
    
    if os.path.exists(pq_path):
        # 读取 parquet
        df = pd.read_parquet(pq_path)
        print(f"读取 {pq_file}, shape={df.shape}")
        
        # 输出 csv
        csv_file = pq_file.replace(".parquet", ".csv")
        csv_path = os.path.join(output_dir, csv_file)
        df.to_csv(csv_path, index=False)
        print(f"已保存为 {csv_path}")
    else:
        print(f"{pq_file} 不存在，跳过")
