# Dataset Overview

**Dataset Source:** Virtual Cell Challenge 2025  
**Input file:** `adata_Training.h5ad`  
**Format:** AnnData (single-cell standard)

**Dataset shape**
- **221,273 cells × 18,080 genes**
- **150 perturbation target genes**
- **38176 non-targeting/control cells**

**Annotations (obs)**
- `target_gene`: perturbed gene  
- `guide_id`: specific gRNA  
- `batch`: sequencing batch

**Gene info (var)**
- `gene_id`: stable ID for each gene

> **Speaker notes:**  
> This study uses the perturbation-based single-cell dataset provided by the Virtual Cell Challenge.  
> It contains ~221k cells, ~18k genes and 150 perturbed (target) genes. Metadata includes gRNA, batch and perturbation labels used as supervision.

---

# Data Distribution

## 1 Severe imbalance in perturbation samples
- Some perturbations have **>10,00 cells**
- Many perturbations have **<500 cells**

## 2 Strong batch effects
- Multiple sequencing batches present
- scRNA-seq batch effects can dominate biological signal

## 3 Large variation in DEG counts across perturbations
- Different perturbations yield drastically different numbers of DEGs  
- Some targets show only a few hundred DEGs (e.g., **MED13: 204 DEGs, 70 up / 134 down**)  
- Others show thousands (e.g., **KDM1A: 3053 DEGs, 2194 up / 859 down**)  
- This creates heterogeneity in signal strength and complicates model training and evaluation

> **Speaker notes:**  
> Beyond sample imbalance and batch effects, perturbations differ dramatically in the strength of their transcriptional response — some produce only modest changes while others induce thousands of DEGs. This wide dynamic range makes it challenging for models to generalize across perturbations with very different effect sizes.

---

# Data Preprocessing

1. **Normalization**  
   Makes gene expression comparable across cells.  

2. **Log1p Transformation**  
   Stabilizes variance and reduces the dominance of highly expressed genes.

3. **Highly Variable Gene (HVG) Selection**  
   Identify the top **2,000 genes** with the highest normalized dispersion.

4. **Output Matrix**  
   Store processed features in `obsm["X_hvg"]` (shape: *n_cells × 2000*).

> **Speaker notes:**  
> We perform standard scRNA-seq preprocessing (normalization + HVG selection), compute virtual-cell embeddings for STATE, and use ESM2 or similar embeddings for perturbations. Splits are by perturbation gene to avoid leakage.

---

# Dataset Split Strategy

## Why split by perturbation gene?
- To evaluate model **generalization to unseen perturbations** (zero-shot)

## Example split (by pert. gene)
- **Training:** 64% of perturbation genes  
- **Validation:** 16%  
- **Test:** 30-gene test set 20%  
- Controls (non-targeting) included in all splits  
- **No overlap** of perturbation genes between splits

## Benefits
- Prevents label leakage across splits
- Tests zero-shot transfer capability
- Ensures fair evaluation of perturbation prediction

> **Speaker notes:**  
> We split at the perturbation-gene level (not by cells). This avoids leakage and lets us measure whether a model can predict effects for genes it never saw during training.

---
# 🧬 STATE Model Reproduction — Summary

---

## 1️⃣ Core Idea of STATE
- **SE (State Embedding):**  
  Learns biologically informed cell embeddings using gene-level ESM2 features + Transformer.
- **ST (State Transition):**  
  Learns how a perturbation transforms control → perturbed cells using a Transformer + **MMD loss**.

**Goal:** Model how genetic perturbations reshape single-cell transcriptomes.

---

## 2️⃣ Model Architecture (High-Level)

### **SE Module**
- Input: top **2048 expressed genes** per cell  
- Gene embeddings from **ESM2 → 672-dim**  
- Add expression-value soft-binning (10 bins)  
- Sequence = genes + **[CLS]** + **[DS]**  
- Transformer encoder → **cell embedding (682-dim)**

### **ST Module**
- Inputs:
  - control cell embedding  
  - perturbation embedding  
  - batch embedding  
- Project to shared hidden dim  
- Transformer backbone models “state shift”  
- Output: predicted perturbed expression  
- Loss: **MMD** between predicted vs real distributions

---
---

## 🌐 总体结构概览

STATE 模型分为两个主要模块：

| 模块 | 名称 | 作用 | 训练方式 |
|------|------|------|-----------|
| 🧠 **SE (State Embedding)** | 学习单细胞的表示 | 自监督训练（self-supervised） | 生成 cell embedding |
| 🔁 **ST (State Transition)** | 学习扰动如何改变细胞状态 | 有监督训练（supervised） | 实现 control → perturbed 转换 |

---

## 🧩 SE 模块（State Embedding, §4.4）

SE 模型用于学习每个细胞的高维 embedding。其核心思路是：
> 利用蛋白语言模型（ESM2）将每个基因转化为语义向量，再通过 Transformer 建模基因间关系，从而获得细胞级别的表征。

### 🔬 结构步骤与论文对应

| 步骤 | 论文出处 | 说明 | 输出维度 |
|------|-----------|------|-----------|
| 1️⃣ **基因层 embedding** | Eq. (23) | 用 **ESM2** 生成每个基因 5120 维 embedding，经线性层压缩至 `h=672` | 每个基因 → 672维 |
| 2️⃣ **细胞层序列** | Eq. (24) | 每个细胞选取 **前 2048 个高表达基因**，加上 `[CLS]` 和 `[DS]` 两个 token，形成序列长度 **2050 (L+2)** | (2050 × 672) |
| 3️⃣ **表达权重嵌入** | Eq. (25)–(27) | 通过 soft binning 将表达值映射至 10 个 bins，得到表达 embedding，加到对应基因 embedding 上 | 不改变维度 (672) |
| 4️⃣ **Transformer 编码** | Eq. (28)–(30) | 用 Transformer 捕捉上下文依赖；取 `[CLS]` 位置输出作为细胞整体表示 | `e_cls ∈ R^672` |
| 5️⃣ **Dataset embedding 拼接** | Eq. (31) | 将 `[DS]` token 的 embedding（降维为 10）拼接到 `e_cls` 后，得到最终细胞表示 | `z_cell ∈ R^(672+10)=R^682` |

---

## 🔁 ST 模块（State Transition, §4.3）

ST 模块学习如何将“control 细胞”转换为“perturbed 细胞”。

### ✳️ 输入组成

| 输入变量 | 含义 | 形状 | 说明 |
|-----------|------|------|------|
| `X_ctrl` | control 细胞表达矩阵 | B × S × G | 控制组细胞的表达（或 embedding） |
| `Z_pert` | perturbation embedding | B × S × D_pert | 扰动条件（如基因敲除） |
| `Z_batch` | batch embedding | B × S × D_batch | 批次信息（技术噪声控制） |

参数定义：
| 符号 | 含义 | 典型值 |
|-------|--------|---------|
| G | 基因数 | 18,080 |
| S | 每个 set 的 cell 数 |  |
| B | batch 数（mini-batch size） |  |
| dh | hidden dim |  |
| D_pert | 扰动 embedding 维度（ESM2 压缩后） |  |

---

### ⚙️ 核心计算流程

1️⃣ **嵌入层**  
所有输入映射到相同 hidden dim：
$$
H = H_{cell} + H_{pert} + H_{batch}
$$
- $ H_{cell} = f_{cell}(X_{ctrl}) $， ( → dh)  
- $ H_{pert} = f_{pert}(Z_{pert}) $， ( → dh)  
- $ H_{batch} = f_{batch}(Z_{batch}) $， ( → dh)

2️⃣ **Transformer Backbone**  
$$
O = H + f_{ST}(H)
$$
- $ f_{ST} $：transformer 层，捕捉 perturbation 对细胞分布的影响。

3️⃣ **输出层**  
线性层将 hidden 表示映射回基因表达空间：  
$$
\hat{X}_{target} = f_{recon}(O) = O W_{recon} + b_{recon}
$$

输出维度：  
$$
\hat{X}_{target} \in \mathbb{R}^{B \times S \times G}
$$

---

### 📉 损失函数：MMD（Maximum Mean Discrepancy）

ST 模型使用 MMD loss 衡量预测分布与真实扰动分布的差异
---

## 3️⃣ Reproduction Workflow

### **Data Preparation**
- Load VCC training data  
- log-transform   

### **ST Module**
 

---


