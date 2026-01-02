## Tahoe-100M Drug Perturbation Experiments

### Overview

- **Tahoe-100M**: largest single-cell perturbation dataset  
- **Cells**: 100,648,790  
- **Samples**: 1,344 (drug × concentration × plate)  
- **Drugs**: 380 chemical compounds  
- **Cell lines**: ~50 cancer cell lines, 14 plates  
- **Perturbation type**: chemical drugs at multiple concentrations + control    

**Cell-level metadata (`obs`):**  
- `plate`: experimental plate ID       
- `phase`: cell cycle phase  
- `S_score`, `G2M_score`: cell cycle scores
- ...... 

---

### Data Extraction

- Original Tahoe lacks **drug target annotations**  
- We obtained **drug-target mappings** from Hugging Face for supervision & evaluation  
- Selected subset for experiments: **52 drugs**  
- Gene matrix: 62,710 genes total; subset has ~31,597 genes  

### Biological Observations

- Target gene expression often **not significantly changed**  
  - Time after drug treatment unknown  
  - Target may affect **protein activity**, not mRNA  
  - Downstream pathway genes often more significant  
- Drugs may have **multiple targets**, behavior varies across cell lines  
- Tahoe reflects **real-world drug target discovery** challenges

---

### Data Split Strategy

- drug level split

- **Single-target drugs**: kept in training & validation  
- **Multi-target drugs**:  
  - Enter test set only if all targets seen in training  
  - Test multi-target drugs have **target-count limit**  
  - High-target drugs kept in training  
- Overlap of target genes between drugs is limited  

---

### Modeling & Evaluation

- Analyses:  
  - Target gene in **top 5/10/50/100 DE genes**  
  - **Logistic Regression baseline** 
 
- Result: baseline cannot detect targets  
- Challenges: complex data, limited sample size, unknown drug effect timing  
- Current approach: different cell lines combined; **future work needed**  

---
