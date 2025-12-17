## 1. 研究动机（Motivation）

我们尝试将**靶点预测（target prediction）**作为扰动预测的下游评价出口。  
通过引入一个**低维、因果意义更明确的任务**，将其作为扰动预测质量的代理指标，从而缓解高维基因表达噪声带来的评价偏差。

---

## 2. 核心思想（Key Idea）

如果一个扰动预测模型能够真实地刻画基因扰动对细胞命运（cell fate）的影响，那么它不仅应当能够预测扰动后的表达变化，还应当支持**反向扰动预测（reverse perturbation prediction）**：

> 在给定扰动后细胞状态的情况下，推断出导致该状态的基因扰动来源。

这一任务在概念上与**药物靶点发现（drug target discovery）**高度相关，因为其目标正是识别驱动细胞状态变化的关键基因。

---

## 3. 相关工作：scGPT

scGPT 在其工作中提出了一种 *in silico* 的反向扰动预测任务，并在 **Norman Perturb-seq** 数据集上进行了验证。  
实验结果表明，scGPT 学到的细胞表示能够在基因组合扰动的复杂空间中，有效地恢复出真实的扰动基因，从而展示了其在靶点推断方面的潜力。

---

## 4. 实验设计（Experimental Design）

### 数据集
- **Norman Perturb-seq** 数据集  
- 选取 20 个基因，共 **210 种扰动组合**

### 任务形式
- 基于扰动结果的 **Top-K 扰动来源检索（Top-K retrieval）**

### 对比方法（Baselines）
- 基于差异表达基因（DEG）的简单启发式方法  
- （可选）GEARS

### 评价指标
- **正确检索（Correct retrieval）**：预测结果与真实扰动完全一致  
- **相关检索（Relevant retrieval）**：预测结果中至少包含真实扰动组合中的一个基因

---

## 5. 与靶点预测的关系（Connection to Target Prediction）

反向扰动预测可以被视为一种**靶点预测的代理任务**。  
它将原本依赖于高维基因表达相似度的比较问题，转化为一个**低维、可验证的判断问题**：

> “该基因是否是导致细胞状态变化的潜在驱动因子？”

---

## 12.14

下载Tahoe-100M metadata
```
mkdir -p tahoe_metadata

gsutil -m cp -r \
  gs://arc-ctc-tahoe100/2025-02-25/metadata \
  ./tahoe_metadata
```

---

## 12.16

Tahoe-100M中有380种药物

ChEMBL:
drug → 已知靶点基因
```
e.g. drug: ADAGRASIB
target_gene: KRAS
mechanism: inhibitor
KRAS G12C 抑制剂
治疗肺癌
```

---

## 12.17

```
['smp_1588' 'smp_1684' 'smp_1780' 'smp_1876' 'smp_1972' 'smp_2068'
 'smp_2152' 'smp_2164' 'smp_2248' 'smp_2260' 'smp_2344' 'smp_2356'
 'smp_2452' 'smp_2548' 'smp_2644' 'smp_2667' 'smp_2668' 'smp_2836']
Number of samples with ADAGRASIB: 18
plate    cell_line
plate8   CVCL_0546    12447
plate13  CVCL_0546    11110
plate2   CVCL_0459     9590
plate13  CVCL_0459     9033
plate2   CVCL_0546     8892
plate9   CVCL_0459     8689
         CVCL_0546     8463
plate8   CVCL_0459     8204
plate12  CVCL_0546     8013
plate13  CVCL_0480     7992
plate8   CVCL_0480     7799
plate7   CVCL_0546     7730
plate6   CVCL_0546     7698
plate7   CVCL_0459     7695
plate10  CVCL_0459     7503
         CVCL_0546     7488
plate8   CVCL_1693     7373
plate11  CVCL_0459     7180
plate4   CVCL_0546     6931
plate7   CVCL_0480     6666
dtype: int64
```
CVCL_0546 结直肠腺癌细胞，KRAS G12V，不是 G12C

CVCL_0459 肺癌细胞

```
mkdir -p tahoe_plate2

gsutil cp \
gs://arc-ctc-tahoe100/2025-02-25/h5ad/plate2_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad \
./tahoe_plate2

```