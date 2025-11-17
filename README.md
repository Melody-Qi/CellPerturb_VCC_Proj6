# VCC Training Data Statistics

文件夹 (`statistics/`) 包含对 VCC 训练数据 (`adata_Training.h5ad`)(`STATE/preprocessed_data/preprocessed_training_data_2000.h5ad`) 的初步数据分析

`STATE/preprocessed_data/preprocessed_training_data_2000.h5ad`由以下指令得到

```
cd ./STATE || exit 1
state tx preprocess_train \
  --adata vcc_data/adata_Training.h5ad \
  --output preprocessed_data/preprocessed_training_data_2000.h5ad \
  --num_hvgs 2000
```
这一步会做：

- Normalize (sc.pp.normalize_total)
- Log1p (Log transform)
- 选取高变基因 (sc.pp.highly_variable_genes) 2000
- 输出 .obsm['X_hvg']，供模型训练使用。
---

## 📌 1️⃣ `training_data_statistics.py`

### 🔍 功能
对全体 **221,273 个细胞 × 18,080 个基因**做统计，包括：

- 每个基因的平均表达 / 最大值 / 最小值  
- 每个基因在多少细胞中有非零表达  
- 总体基因表达分布图    

### 📈 输出示例
```
读取数据：../STATE/vcc_data/adata_Training.h5ad
AnnData object with n_obs × n_vars = 221273 × 18080
obs: target_gene, guide_id, batch
var: gene_id

平均表达值均值： 3.1473236
平均最大表达值： 32.037113
平均最小值表达值： 0.007522124
基因完全不表达的比例： 0.0002212389
```
```
读取数据：../STATE/preprocessed_data/preprocessed_training_data_2000.h5ad
AnnData object with n_obs × n_vars = 221273 × 18080
    obs: 'target_gene', 'guide_id', 'batch'
    var: 'gene_id', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'hvg', 'log1p'
    obsm: 'X_hvg'
obs keys: Index(['target_gene', 'guide_id', 'batch'], dtype='object')
var keys: Index(['gene_id', 'highly_variable', 'means', 'dispersions',
       'dispersions_norm'],
      dtype='object')
平均表达值均值： 0.7703777
平均最大表达值： 2.6602652
平均最小值表达值： 0.005208211
基因完全不表达的比例： 0.00022123893805309734
保存统计文件: preprocessed_training_data_2000/preprocessed_training_data_2000_gene_statistics.csv
保存统计图: preprocessed_training_data_2000/plots/preprocessed_training_data_2000_gene_mean_distribution.png
```
| 字段                 | 含义                 | 作用               |
| ------------------ | ------------------ | ---------------- |
| `means`            | 每个基因在全部细胞中的平均表达    | 用于评估基因表达水平       |
| `dispersions`      | 原始变异度（基因表达方差 / 均值） | 衡量基因表达变化大小       |
| `dispersions_norm` | 标准化后的变异度（去除表达量偏置）  | 用于排序并筛选 HVGs     |
| `highly_variable`  | True / False       | 是否被标记为高变基因 (HVG) |

| 情况                          | 含义                          |
| --------------------------- | --------------------------- |
| 高 mean & 高 dispersions_norm | 信息量高的重点基因                   |
| 高 mean & 低 dispersions_norm | 高表达但不变 = housekeeping genes |
| 低 mean & 高 dispersions_norm | 稀有信号（可能有用，也可能是噪声）           |
| 低 mean & 低 dispersions_norm | 基本无信息基因                     |


### 🧪 生成的 CSV
`gene_statistics.csv`
`preprocessed_training_data_2000_gene_statistics`（前几行示例）：
| gene     | mean | max  | min  | nonzero_fraction |
|----------|------|------|------|-----------------|
| SAMD11   |  |  |   |             |
| NOC2L    |  |  |   |             |
| ...      | ...  | ...  | ...  | ...             |

### 📊 输出图示 (plots/)
- **18080 基因平均表达分布** 
![illustration](./statistics/plots/gene_mean_distribution.png) 
  说明：大多数基因表达低，少数高表达（符合 scRNA-seq 偏态分布）

---

## 📌 2️⃣ `gene_delta_statistics.py`

### 📈 数据概况
- 总细胞: 221,273  
- control: 38,176  
- perturbed: 183,097  
- Control 中永不表达的基因数量：239  

### 🧪 生成的 CSV
`gene_delta_all_sorted.csv`（按 `abs_delta` 从大到小排序）：
| gene   | ctrl_mean | pert_mean | delta  | abs_delta |
|--------|-----------|-----------|--------|-----------|
| CD24   |      |       |   |       |
| HSPA8  |       |       |   |       |
| ...    | ...       | ...       | ...    | ...       |

- `gene_delta_top50.csv`：Top-K 差异基因，直接用于分析或模型 debug  
- `gene_delta_with_nonzero_fraction.csv`：附加 control 非零表达比例，有助于标记永不表达基因  

### 📊 输出图示 (plots/)
- **基因表达变化分布（delta = pert_mean - ctrl_mean）**  
![illustration](./statistics/plots/gene_delta_distribution.png)
  说明：扰动信号极稀疏 → 大部分基因变化 ≈ 0，少量显著变化  

- **Control vs Perturbed 散点图（log scale）**  
![illustration](./statistics/plots/ctrl_vs_pert_scatter.png)
  说明：点集中在对角线 → 大多数基因不受扰动影响  

- **Top 50 差异基因条形图**  
![illustration](./statistics/plots/top50_delta_barh.png)
  说明：模型应该重点学习这些 gene

---

## 📎 脚本使用方法

**全基因统计：**
```bash
cd statistics
python training_data_statistics.py -i ../STATE/vcc_data/adata_Training.h5ad
