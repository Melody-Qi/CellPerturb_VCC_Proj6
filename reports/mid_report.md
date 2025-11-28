# Dataset Overview

**Dataset Source:** Virtual Cell Challenge 2025  
**Input file:** `adata_Training.h5ad`  
**Format:** AnnData (single-cell standard)

**Dataset shape**
- **221,273 cells × 18,080 genes**
- **150 perturbation target genes**

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

# Slide 2 — Data Distribution: Why This Task Is Hard

## 1 Severe imbalance in perturbation samples
- Some perturbations have **>10,000 cells**
- Many perturbations have **<300 cells**
- → **Long-tail distribution**: hard for models to learn rare perturbations

*(Recommend plotting: bar plot of `target_gene` counts; show top 20 + "others")*

## 2 Strong batch effects
- Multiple sequencing batches present
- scRNA-seq batch effects can dominate biological signal
- → Model must learn **batch-invariant** representations

*(Recommend plotting: histogram/pie chart of `batch` counts)*

## 3 Control vs perturbed imbalance
- Control (`non-targeting`) cells are limited
- Perturbed cells dominate but are unevenly distributed
- → Mean-based evaluation is sensitive to this imbalance

*(Recommend reporting counts: #control cells vs #perturbed cells)*

> **Speaker notes:**  
> The dataset poses three core challenges: severe class imbalance across perturbations, clear batch effects, and limited control cells — all of which make prediction and fair evaluation difficult.

---

# Slide 3 — Data Preprocessing

## 1 Normalization
- Apply `log1p` or CPM/TMM normalization
- Rationale: scRNA-seq counts are highly skewed; normalization stabilizes training

## 2 Highly Variable Gene (HVG) selection
- Typical choice: **2k–5k HVGs**
- Purpose: remove noisy genes, reduce dimensionality and compute cost
- Use HVGs as `X_hvg` input to models

## 3 Virtual cell embedding (STATE-specific)
- Use VAE / cell encoder to map per-cell expression → low-dim latent
- These embeddings form the model input (virtual cells / cell sets)

## 4 Perturbation feature processing
- Encode `target_gene` via:
  - ESM2 perturbation embeddings (provided by challenge), OR
  - one-hot vectors, OR
  - learned embeddings
- Normalize or standardize perturbation features as needed

## 5 Data split principle
- **Split by perturbation gene** (not by cell)
- Prevents leakage and allows evaluation of zero-shot generalization

> **Speaker notes:**  
> We perform standard scRNA-seq preprocessing (normalization + HVG selection), compute virtual-cell embeddings for STATE, and use ESM2 or similar embeddings for perturbations. Splits are by perturbation gene to avoid leakage.

---

# Slide 4 — Dataset Split Strategy

## Why split by perturbation gene?
- To evaluate model **generalization to unseen perturbations** (zero-shot)

## Example split (by pert. gene)
- **Training:** 80% of perturbation genes  
- **Validation:** 10%  
- **Test:** 10%  
- Controls (non-targeting) included in all splits  
- **No overlap** of perturbation genes between splits

## Benefits
- Prevents label leakage across splits
- Tests zero-shot transfer capability
- Ensures fair evaluation of perturbation prediction

> **Speaker notes:**  
> We split at the perturbation-gene level (not by cells). This avoids leakage and lets us measure whether a model can predict effects for genes it never saw during training.

---

# Appendix — Suggested Figures and Commands

## Suggested figures to create (quick summary)
1. Bar plot: `target_gene` counts (sorted, top-20 + others)
2. Histogram / pie chart: `batch` distribution
3. Table: counts of control vs perturbed cells
4. HVG selection: scree plot or histogram of gene variance; show chosen HVG count

## Quick preprocessing CLI examples
```bash
# select 2000 HVGs (STATE preprocessing example)
state tx preprocess_train \
  --adata ./vcc_data/adata_Training.h5ad \
  --output ./preprocessed_data/preprocessed_training_data_2000.h5ad \
  --num_hvgs 2000
