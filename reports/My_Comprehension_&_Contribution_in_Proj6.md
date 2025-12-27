# Predicting the response of cellular transcriptome to gene perturbations

## Completing the VCC Competition

SOTA: scGPT **STATE**

### preprocess_train

- Normalization (sc.pp.normalize_total)

- Log1p transformation (log transform)

- Selection of highly variable genes (HVGs) using sc.pp.highly_variable_genes (top 2000)

- Outputs .obsm['X_hvg'] for model training
---

### train

The data and configuration files used for training are as follows:

| Item                  | File path                                        | Description                                     |
| --------------------- | ------------------------------------------------ | ----------------------------------------------- |
| Training data         | `competition_support_set/competition_train.h5ad` | Training set (including control + perturbation) |
| Config file           | `competition_support_set/starter.toml`           | Specifies data loading and column mappings      |
| Perturbation features | `competition_support_set/ESM2_pert_features.pt`  | ESM2 embeddings for gene perturbations          |
| Output directory      | `./checkpoints/first_run`                        | Stores model checkpoints and logs               |
---

### Reading the STATE Paper

(2025). Single-cell Transformer for Adaptive Transcriptomic Effects (STATE)
Relevant sections: Methods §4.3 (ST) & §4.4 (SE)

---

The STATE model consists of two main modules:

🧩SE Module (State Embedding, §4.4)

The SE module learns high-dimensional embeddings for individual cells.

The core idea is to use a protein language model (ESM2) to encode genes into semantic vectors and apply a Transformer to model gene–gene relationships.

| Step                 | Description                                          |
| -------------------- | ---------------------------------------------------- |
| Gene embedding       | ESM2 generates 5120-dim embeddings, projected to 672 |
| Cell sequence        | Top 2048 expressed genes + `[CLS]`, `[DS]`           |
| Expression embedding | Soft binning into 10 bins                            |
| Transformer          | `[CLS]` output represents the cell                   |
| Dataset embedding    | `[DS]` concatenated → final 682-dim embedding        |

---

🔁 ST Module (State Transition, §4.3)

The ST module learns how to transform a control cell state into a perturbed cell state.

✳️ Inputs

| Input     | Meaning                 | Shape           | Description                                    |
| --------- | ----------------------- | --------------- | ---------------------------------------------- |
| `X_ctrl`  | Control cell expression | B × S × G       | Expression (or embeddings) of control cells    |
| `Z_pert`  | Perturbation embedding  | B × S × D_pert  | Perturbation condition (e.g., gene knockout)   |
| `Z_batch` | Batch embedding         | B × S × D_batch | Batch information (to control technical noise) |

Parameter definitions:
| Symbol | Meaning                                                   | Typical value  |
| ------ | --------------------------------------------------------- | -------------- |
| G      | Number of genes                                           | 18,080         |
| S      | Number of cells per set                                   | —              |
| B      | Batch size (mini-batch)                                   | Depends on GPU |
| dh     | Hidden dimension                                          | —              |
| D_pert | Perturbation embedding dimension (after ESM2 compression) | —              |

---

⚙️ Core Computation Flow

All inputs are projected into the same hidden dimension:
$$
H = H_{cell} + H_{pert} + H_{batch}
$$
- $ H_{cell} = f_{cell}(X_{ctrl}) $， ( → dh)  
- $ H_{pert} = f_{pert}(Z_{pert}) $， ( → dh)  
- $ H_{batch} = f_{batch}(Z_{batch}) $， ( → dh)

**Transformer Backbone**  
$$
O = H + f_{ST}(H)
$$
- $ f_{ST} $：Transformer layers modeling how perturbations affect cell distributions.

**Output layer**  
A linear layer maps hidden representations back to gene expression space:  
$$
\hat{X}_{target} = f_{recon}(O) = O W_{recon} + b_{recon}
$$

Output shape：  
$$
\hat{X}_{target} \in \mathbb{R}^{B \times S \times G}
$$

---

📉 Loss Function：MMD（Maximum Mean Discrepancy）

The ST module uses MMD loss to measure the discrepancy between the predicted and true perturbed distributions:

$$
L_{MMD}(\hat{X}_{target}, X_{target}) = \frac{1}{B} \sum_{b=1}^{B} MMD^2(\hat{X}_b, X_b)
$$

where：
$$
MMD^2(\hat{X}, X) = \frac{1}{S^2} \sum_{i,j} [k(\hat{x}_i,\hat{x}_j) + k(x_i,x_j) - 2k(\hat{x}_i, x_j)]
$$
The kernel used is the **energy distance kernel**：
$$
k(u,v) = -\|u - v\|_2
$$

This loss encourages the generated perturbation distribution to match the statistical properties of the real distribution.

---

### statistics

Statistical analysis over **221,273 cells × 18,080 genes**, including:

- Mean / max / min expression per gene

- Number of cells with non-zero expression per gene

- Global gene expression distribution


🧪 Generated CSV

example rows:
| gene     | mean | max  | min  | nonzero_fraction |
|----------|------|------|------|-----------------|
| SAMD11   |  |  |   |             |
| NOC2L    |  |  |   |             |
| ...      | ...  | ...  | ...  | ...             |

📈 Output examples
```
Loading data: ../STATE/vcc_data/adata_Training.h5ad
AnnData object with n_obs × n_vars = 221273 × 18080
obs: target_gene, guide_id, batch
var: gene_id

Mean of gene means: 3.1473236
Mean of gene maxima: 32.037113
Mean of gene minima: 0.007522124
Fraction of genes never expressed: 0.0002212389
```
```
Loading data: ../STATE/preprocessed_data/preprocessed_training_data_2000.h5ad
AnnData object with n_obs × n_vars = 221273 × 18080
    obs: 'target_gene', 'guide_id', 'batch'
    var: 'gene_id', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'hvg', 'log1p'
    obsm: 'X_hvg'

obs keys: Index(['target_gene', 'guide_id', 'batch'], dtype='object')
var keys: Index(['gene_id', 'highly_variable', 'means', 'dispersions',
       'dispersions_norm'], dtype='object')

Mean of gene means: 0.7703777
Mean of gene maxima: 2.6602652
Mean of gene minima: 0.005208211
Fraction of genes never expressed: 0.00022123893805309734

Saved statistics file:
preprocessed_training_data_2000/preprocessed_training_data_2000_gene_statistics.csv

Saved plot:
preprocessed_training_data_2000/plots/preprocessed_training_data_2000_gene_mean_distribution.png
```
Field descriptions
| Field              | Meaning                                         | Usage                                  |
| ------------------ | ----------------------------------------------- | -------------------------------------- |
| `means`            | Mean expression per gene across all cells       | Assess gene expression level           |
| `dispersions`      | Raw dispersion (variance / mean)                | Measure expression variability         |
| `dispersions_norm` | Normalized dispersion (expression bias removed) | Used for HVG ranking                   |
| `highly_variable`  | True / False                                    | Whether the gene is selected as an HVG |

Interpretation
| Case                              | Meaning                                          |
| --------------------------------- | ------------------------------------------------ |
| High mean & high dispersions_norm | Highly informative genes                         |
| High mean & low dispersions_norm  | Highly expressed but stable (housekeeping genes) |
| Low mean & high dispersions_norm  | Rare signals (potentially useful or noisy)       |
| Low mean & low dispersions_norm   | Largely uninformative genes                      |


📊 (plots/)
- **Mean gene expression distribution across 18,080 genes** 
![illustration](../statistics/plots/gene_mean_distribution.png) 
  Explanation: Most genes have low expression, with a small number of highly expressed genes, consistent with the skewed distribution of scRNA-seq data.

---

📌`gene_delta_statistics.py` 

📊 (plots/)
- **Gene expression change distribution (delta = pert_mean - ctrl_mean)**  
![illustration](../statistics/plots/gene_delta_distribution.png)
  Explanation: Perturbation effects are highly sparse — most genes show changes close to zero, with only a small number exhibiting significant shifts.

---

## Dataset Split Strategy

**Goal:** Evaluate generalization to **unseen perturbations (zero-shot)**

- **Training:** 64% of perturbation genes  
- **Validation:** 16% of perturbation genes  
- **Test:** 20% (30-gene zero-shot set)

**Notes:**
- Non-targeting **controls included in all splits**
- **No overlap** of perturbation genes between splits

## Evaluation Results & Reflections

| Model                     | DES ↑    | PDS ↑    | MAE ↓    | Overall ↑ |
|---------------------------|----------|----------|----------|-----------|
| Cell-mean Baseline        | 0.0442   | 0.5167   | 0.1258   | 0.00      |
| Ridge Regression          | 0.1055   | 0.5156   | 0.0895   | 11.76     |
| Random Forest Regression  | 0.0975   | 0.5567   | 0.1007   | 11.27     |
| State Model               | 0.0837   | 0.4478   | 0.1473   | 1.38      |
| State Model (More Steps)  | 0.0633   | 0.5122   | 0.1050   | 6.18      |
| scGPT (Zero-shot)         | 0.1767   | 0.4833   | 1.4721   | 4.62      |
| scGPT (Finetuned)         | 0.1870   | 0.4833   | 1.2272   | 4.98      |


# 从扰动预测走向下游应用：药物靶点预测