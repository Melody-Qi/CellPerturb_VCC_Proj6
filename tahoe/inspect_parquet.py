import pandas as pd
import pyarrow.parquet as pq
import os
import numpy as np

# 查看生成的数据目录结构
OUT_DIR = "./tahoe_under_10gb_subset"
print("输出目录内容:")
for file in os.listdir(OUT_DIR):
    if file.endswith('.parquet'):
        file_path = os.path.join(OUT_DIR, file)
        file_size = os.path.getsize(file_path) / (1024**2)  # MB
        print(f"  {file} ({file_size:.2f} MB)")

# 读取第一个parquet文件查看数据
parquet_files = [f for f in os.listdir(OUT_DIR) if f.endswith('.parquet')]
if parquet_files:
    first_file = os.path.join(OUT_DIR, parquet_files[0])
    print(f"\n读取文件: {first_file}")
    
    # 读取数据
    df = pq.read_table(first_file).to_pandas()
    
    print(f"数据形状: {df.shape}")
    print(f"\n列名: {list(df.columns)}")
    
    print("\n前5行数据:")
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 50)
    print(df.head())
    
    print("\n数据基本信息:")
    print(df.info())
    
    print("\n数值列统计:")
    print(df.describe())
    
    # 查看y值的分布（类别分布）
    print(f"\n类别分布 (y值):")
    y_counts = df['y'].value_counts().sort_index()
    print(y_counts)
    
    # 查看target_gene的分布
    print(f"\n靶基因分布:")
    target_counts = df['target_gene'].value_counts().head(10)
    print(target_counts)
    
    # 查看drug的分布
    print(f"\n药物分布 (前10):")
    drug_counts = df['drug'].value_counts().head(10)
    print(drug_counts)
    
    # 查看cell_line_id的分布
    print(f"\n细胞系分布 (前10):")
    cell_line_counts = df['cell_line_id'].value_counts().head(10)
    print(cell_line_counts)
    
    # 查看gene_ids和counts的结构
    print(f"\n基因数据示例:")
    for i in range(min(3, len(df))):
        print(f"\n样本 {i}:")
        print(f"  y值: {df.iloc[i]['y']}")
        print(f"  靶基因: {df.iloc[i]['target_gene']}")
        print(f"  药物: {df.iloc[i]['drug']}")
        print(f"  基因数量: {len(df.iloc[i]['gene_ids'])}")
        print(f"  前5个基因ID: {df.iloc[i]['gene_ids'][:5]}")
        print(f"  前5个表达量: {df.iloc[i]['counts'][:5]}")
        
        # 检查基因表达量的统计信息
        counts = np.array(df.iloc[i]['counts'])
        print(f"  表达量统计 - 最小值: {counts.min():.2f}, 最大值: {counts.max():.2f}, 平均值: {counts.mean():.2f}")
    
    # 检查数据完整性
    print(f"\n数据完整性检查:")
    print(f"空值统计:")
    print(df.isnull().sum())
    
    # 检查gene_ids和counts的长度是否一致
    lengths_consistent = all(len(row['gene_ids']) == len(row['counts']) for _, row in df.iterrows())
    print(f"基因ID和表达量长度一致: {lengths_consistent}")
    
    # 检查基因数量的分布
    gene_counts = [len(row['gene_ids']) for _, row in df.iterrows()]
    print(f"每个样本的基因数量 - 最小值: {min(gene_counts)}, 最大值: {max(gene_counts)}, 平均值: {np.mean(gene_counts):.1f}")

else:
    print("没有找到parquet文件")