import torch
from unipert import UniPert

import pandas as pd

unipert = UniPert()

# 1. 读取基因列表
gene_file = "../../vcc_data/gene_names.csv"
with open(gene_file, "r") as f:
    genes = [line.strip() for line in f.readlines()]

# 2. 生成 embedding
# encode_genes 返回 tensor (num_genes, d)
out_embs, invalid_inputs = unipert.enc_gene_ptbgs_from_gene_names(gene_names=genes)

# 3. 保存到 .pt 文件
torch.save(out_embs, "my_unipert_features.pt")

print("✅ 已生成 UniPert embedding 并保存为 my_unipert_features.pt")

# 4. 保存到 CSV 文件，查看 embedding
# 转成 numpy
emb_np = out_embs.numpy()  # shape = (num_genes, d)

# 创建 DataFrame
df = pd.DataFrame(emb_np, index=genes)
df.index.name = "gene_id"

# 保存 CSV
df.to_csv("my_unipert_features.csv")
print("✅ 已生成 CSV 文件 my_unipert_features.csv, shape =", emb_np.shape)
