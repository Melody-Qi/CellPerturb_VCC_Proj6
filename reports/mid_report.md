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

## 3️⃣ Reproduction Workflow

### **Data Preparation**
- Load VCC training data  
- log-transform   

### **ST Module**
 

---

## 4️⃣ Fine-tuning Pipeline

---

