# python my_cell_eval3.py
import json
import warnings
from typing import Optional, Dict, Union
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import anndata as ad
from hpdex import parallel_differential_expression

def _get_matrix_row_mean(adata, mask):
    """
    返回 mask（布尔数组)对应细胞的每个基因的平均表达（1D numpy float array)。
    兼容稀疏矩阵。
    """
    X = adata.X
    if sp.issparse(X):
        sub = X[mask]
        # mean(axis=0) 返回 1xG sparse/dense matrix
        mean_vec = np.array(sub.mean(axis=0)).ravel()
    else:
        sub = X[mask, :]
        mean_vec = np.asarray(sub.mean(axis=0)).ravel()
    return mean_vec.astype(float)

def _extract_expr_vectors_for_gene(adata, mask):
    """返回 mask 对应细胞在所有基因上的表达子矩阵（稀疏或密集)"""
    X = adata.X
    if sp.issparse(X):
        sub = X[mask].toarray()
    else:
        sub = X[mask, :]
    return sub  # shape (n_cells, n_genes)

def _wilcoxon_de_genes(adata, group_mask, ntc_mask, gene_names, fdr=0.05, sample_n:int=1000):
    """
    对每个基因做 Wilcoxon rank-sum（Mann-Whitney U)检验：group vs ntc。
    为控制时间，当每组细胞非常多时，会进行随机抽样（sample_n)。
    返回：DataFrame 包含 index=gene_names, columns=['pval','logfc','mean_group','mean_ntc']
    """
    n_group = int(group_mask.sum())
    n_ntc = int(ntc_mask.sum())
    if n_group < 2 or n_ntc < 2:
        # 无法检验
        return pd.DataFrame(index=gene_names, data={
            'pval': np.ones(len(gene_names)),
            'logfc': np.zeros(len(gene_names)),
            'mean_group': np.zeros(len(gene_names)),
            'mean_ntc': np.zeros(len(gene_names))
        })

    # 抽样索引
    rng = np.random.default_rng(0)
    if sample_n is not None:
        if n_group > sample_n:
            g_idx = rng.choice(np.nonzero(group_mask)[0], sample_n, replace=False)
            group_mask_sample = np.zeros_like(group_mask, dtype=bool)
            group_mask_sample[g_idx] = True
        else:
            group_mask_sample = group_mask
        if n_ntc > sample_n:
            n_idx = rng.choice(np.nonzero(ntc_mask)[0], sample_n, replace=False)
            ntc_mask_sample = np.zeros_like(ntc_mask, dtype=bool)
            ntc_mask_sample[n_idx] = True
        else:
            ntc_mask_sample = ntc_mask
    else:
        group_mask_sample = group_mask
        ntc_mask_sample = ntc_mask

    # 取表达矩阵
    group_mat = _extract_expr_vectors_for_gene(adata, group_mask_sample)  # (n_g, G)
    ntc_mat = _extract_expr_vectors_for_gene(adata, ntc_mask_sample)      # (n_n, G)

    G = group_mat.shape[1]
    pvals = np.ones(G)
    mean_group = group_mat.mean(axis=0)
    mean_ntc = ntc_mat.mean(axis=0)
    # logFC 使用 log2(mean_group + 1) - log2(mean_ntc + 1)
    logfc = np.log2(mean_group + 1) - np.log2(mean_ntc + 1)

    # 对每个基因做 Mann-Whitney U（两侧)
    for gi in range(G):
        try:
            u = mannwhitneyu(group_mat[:, gi], ntc_mat[:, gi], alternative='two-sided')
            pvals[gi] = u.pvalue if u.pvalue is not None else 1.0
        except Exception:
            pvals[gi] = 1.0

    # FDR 校正
    reject, pvals_adj, _, _ = multipletests(pvals, alpha=fdr, method='fdr_bh')
    df = pd.DataFrame({
        'pval': pvals,
        'pval_adj': pvals_adj,
        'significant': reject,
        'logfc': logfc,
        'mean_group': mean_group,
        'mean_ntc': mean_ntc
    }, index=gene_names)
    return df


def _load_de_table(df_or_path: Optional[Union[str, pd.DataFrame]]):
    """
    读取或标准化一个 DE 结果表（来自 CSV 或 DataFrame）。
    返回 standardized DataFrame with columns: ['target','feature','fdr','log2fc']
    若参数为 None，返回 None。
    """
    if df_or_path is None:
        return None
    # 如果是文件路径，读取
    if isinstance(df_or_path, str):
        df = pd.read_csv(df_or_path)
    elif isinstance(df_or_path, pd.DataFrame):
        df = df_or_path.copy()
    else:
        raise ValueError("true_de_df / pred_de_df must be a file path or a pandas.DataFrame (or None).")

    # 统一小写列名便于匹配
    colmap = {c.lower(): c for c in df.columns}
    # identify columns for fdr and log2fc
    fdr_col = None
    for cand in ['fdr', 'pval_adj', 'padj', 'adj_pval', 'adj_pvalue']:
        if cand in colmap:
            fdr_col = colmap[cand]
            break
    if fdr_col is None:
        # fallback: use p_value and then will BH adjust later (less ideal but still)
        if 'p_value' in colmap or 'pval' in colmap or 'p_value' in [c.lower() for c in df.columns]:
            # leave fdr_col None -> caller might choose to compute BH, but here we set fdr to pval to be safe
            pass

    logfc_col = None
    for cand in ['log2_fold_change', 'log2fc', 'logfc', 'log2_foldchange']:
        if cand in colmap:
            logfc_col = colmap[cand]
            break

    target_col = None
    for cand in ['target', 'target_gene', 'perturbation', 'group']:
        if cand in colmap:
            target_col = colmap[cand]
            break
    if target_col is None:
        raise ValueError("DE table has no target column (expected 'target'/'target_gene' etc.).")

    feature_col = None
    for cand in ['feature', 'gene', 'gene_name', 'geneid']:
        if cand in colmap:
            feature_col = colmap[cand]
            break
    if feature_col is None:
        raise ValueError("DE table has no feature/gene column (expected 'feature'/'gene' etc.).")

    # build standardized df
    std = pd.DataFrame()
    std['target'] = df[target_col].astype(str)
    std['feature'] = df[feature_col].astype(str)

    # fdr: if available use it; else try to use pval and BH adjust
    if fdr_col is not None:
        std['fdr'] = pd.to_numeric(df[fdr_col], errors='coerce').fillna(1.0)
    else:
        # try pval column
        pcol = None
        for cand in ['p_value', 'pval', 'p.value', 'pvalue', 'p_value']:
            if cand in colmap:
                pcol = colmap[cand]
                break
        if pcol is not None:
            pvals = pd.to_numeric(df[pcol], errors='coerce').fillna(1.0).values
            # perform BH per entire table (conservative but acceptable); note we don't know grouping yet
            _, p_adj, _, _ = multipletests(pvals, method='fdr_bh')
            std['fdr'] = p_adj
        else:
            std['fdr'] = 1.0  # fallback: no significance

    # log2fc: if absent, try to compute from fold_change column if present
    if logfc_col is not None:
        std['log2fc'] = pd.to_numeric(df[logfc_col], errors='coerce').fillna(0.0)
    else:
        # try fold_change
        fc_col = None
        for cand in ['fold_change', 'foldchange', 'fc']:
            if cand in colmap:
                fc_col = colmap[cand]
                break
        if fc_col is not None:
            fc_vals = pd.to_numeric(df[fc_col], errors='coerce').fillna(1.0)
            # avoid log of nonpositive
            with np.errstate(divide='ignore', invalid='ignore'):
                log2fc = np.log2(fc_vals.replace(0, np.nan)).fillna(0.0)
            std['log2fc'] = log2fc
        else:
            std['log2fc'] = 0.0

    return std[['target', 'feature', 'fdr', 'log2fc']]

def compute_vcc_scores(true_adata,
                       pred_adata,
                       baseline_scores: Optional[Dict[str,float]] = None,                       
                       groupby: str = 'target_gene',
                       control_label: str = 'non-targeting',
                       fdr: float = 0.05,
                       sample_n: int = 1000,
                       min_cells_per_group: int = 3,
                       true_de_df: Optional[Union[str,pd.DataFrame]] = None,
                       pred_de_df: Optional[Union[str,pd.DataFrame]] = None):
    """
    计算 DES, PDS, MAE 以及最终 S。
    新增参数：
        true_de_df / pred_de_df: 可选路径或 DataFrame，表示 Wilcoxon 的结果 CSV（columns include target, feature, fdr, log2_fold_change 等）。
                                 若提供，则 DES 的 G_true / G_pred 直接从这些表中读取（使用 fdr阈值），而不再对 AnnData 做 Wilcoxon 计算。
    其余行为与原函数一致（PDS/MAE 没变）。
    """
    # --- 对齐基因（var_names) ---
    genes_true = list(true_adata.var_names)
    genes_pred = list(pred_adata.var_names)
    if genes_true != genes_pred:
        # 取交集并重新索引
        common = [g for g in genes_true if g in genes_pred]
        if len(common) == 0:
            raise ValueError("No common genes between true_adata and pred_adata.")
        true_adata = true_adata[:, common]
        pred_adata = pred_adata[:, common]

    gene_names = list(true_adata.var_names)
    G = len(gene_names)

    # --- DE CSV (若有) 读取并标准化 ---
    true_de_table = _load_de_table(true_de_df)  # may be None
    pred_de_table = _load_de_table(pred_de_df)  # may be None

    # 如果提供了 DE 表，但其中基因/perturbation名与 adata 不匹配，发出警告并用交集
    if true_de_table is not None:
        # only keep features that are in gene_names
        before = true_de_table.shape[0]
        true_de_table = true_de_table[true_de_table['feature'].isin(gene_names)].copy()
        after = true_de_table.shape[0]
        if after < before:
            warnings.warn("Some genes in provided true_de_df not found in AnnData var_names; they were dropped.")
    if pred_de_table is not None:
        before = pred_de_table.shape[0]
        pred_de_table = pred_de_table[pred_de_table['feature'].isin(gene_names)].copy()
        after = pred_de_table.shape[0]
        if after < before:
            warnings.warn("Some genes in provided pred_de_df not found in AnnData var_names; they were dropped.")

    # --- 确定 perturbation 列与分组 ---
    true_groups = true_adata.obs[groupby].astype(str) # groupby = target_gene
    pred_groups = pred_adata.obs[groupby].astype(str)

    # 找到控制(ntc) mask
    ntc_mask_true = (true_groups == control_label).values
    ntc_mask_pred = (pred_groups == control_label).values

    if ntc_mask_true.sum() < min_cells_per_group or ntc_mask_pred.sum() < min_cells_per_group:
        warnings.warn("Control (ntc) cell count is small.")

    # 获取 perturbation 列表（排除 control）
    perturbations = sorted([g for g in true_groups.unique() if g != control_label])
    N = len(perturbations)
    if N == 0:
        raise ValueError("No perturbations found (after excluding control_label).")

    # --- 1) 计算 DES per perturbation ---
    DES_list = []
    for p in perturbations:
        # masks
        g_mask_true = (true_groups == p).values
        g_mask_pred = (pred_groups == p).values

        if g_mask_true.sum() < min_cells_per_group or ntc_mask_true.sum() < min_cells_per_group:
            DES_list.append(np.nan)
            continue

        # 1a) G_true: 如果提供了 true_de_table，则直接从表中取 G_true；否则运行 Wilcoxon
        if true_de_table is not None:
            sub = true_de_table[true_de_table['target'] == p]
            G_true = set(sub.loc[sub['fdr'] <= fdr, 'feature'].tolist())
        else:
            df_true = _wilcoxon_de_genes(true_adata, g_mask_true, ntc_mask_true, gene_names, fdr=fdr, sample_n=sample_n)
            G_true = set(df_true.index[df_true['significant']].tolist())

        # 1b) G_pred: 如果提供了 pred_de_table，则直接从表中取 G_pred；否则运行 Wilcoxon on pred_adata
        if pred_de_table is not None:
            subp = pred_de_table[pred_de_table['target'] == p]
            G_pred = set(subp.loc[subp['fdr'] <= fdr, 'feature'].tolist())
        else:
            if g_mask_pred.sum() < min_cells_per_group or ntc_mask_pred.sum() < min_cells_per_group:
                df_pred = pd.DataFrame(index=gene_names, data={
                    'pval': np.ones(G),
                    'pval_adj': np.ones(G),
                    'significant': np.zeros(G, dtype=bool),
                    'logfc': np.zeros(G),
                    'mean_group': np.zeros(G),
                    'mean_ntc': np.zeros(G)
                })
            else:
                df_pred = _wilcoxon_de_genes(pred_adata, g_mask_pred, ntc_mask_pred, gene_names, fdr=fdr, sample_n=sample_n)
            G_pred = set(df_pred.index[df_pred['significant']].tolist())

        # 若 |G_true| == 0：按原逻辑跳过（记 NaN）
        if len(G_true) == 0:
            DES_list.append(np.nan)
            continue

        # Case 分支：|G_pred| <= |G_true| 或 > 
        if len(G_pred) <= len(G_true):
            inter = len(G_pred & G_true)
            DES_k = inter / len(G_true)
            DES_list.append(DES_k)
        else:
            # > 情况：需要按预测的绝对 log2FC 排序取 top |G_true|
            # 如果 pred_de_table 提供了 log2fc，则使用它；若没有，则尽量使用 df_pred 中计算结果
            if pred_de_table is not None:
                subp_all = pred_de_table[pred_de_table['target'] == p].set_index('feature')
                # 若某些基因在该 perturbation 的表中缺失，补 0
                # 这里确保索引覆盖 gene_names
                # 取 abs(log2fc)
                subp_all = subp_all.reindex(gene_names).fillna({'log2fc': 0.0})
                subp_all['abs_log2fc'] = subp_all['log2fc'].abs()
                topk = int(len(G_true))
                top_genes = set(subp_all.sort_values('abs_log2fc', ascending=False).head(topk).index.tolist())
            else:
                # df_pred 已在上面构造
                df_pred_sorted = df_pred.reindex(gene_names).copy()
                df_pred_sorted['abs_logfc'] = np.abs(df_pred_sorted['logfc'])
                topk = int(len(G_true))
                top_genes = set(df_pred_sorted.sort_values('abs_logfc', ascending=False).head(topk).index.tolist())

            inter = len(top_genes & G_true)
            DES_k = inter / len(G_true)
            DES_list.append(DES_k)

    DES_array = np.array(DES_list, dtype=float)
    DES_mean = np.nanmean(DES_array)  # 忽略 NaN

    # --- 2) 计算 PDS --- (Modifying PDS section)
    # 计算所有 perturbation 的 pseudobulk mean expression (true & pred), 包括 ntc
    true_pseudobulk = {}
    pred_pseudobulk = {}
    for grp in list(true_groups.unique()):
        mask = (true_groups == grp).values
        true_pseudobulk[grp] = _get_matrix_row_mean(true_adata, mask)
    for grp in list(pred_groups.unique()):
        mask = (pred_groups == grp).values
        pred_pseudobulk[grp] = _get_matrix_row_mean(pred_adata, mask)

    for p in perturbations:
        if p not in pred_pseudobulk:
            warnings.warn(f"Perturbation {p} missing in predicted data; using zeros vector.")
            pred_pseudobulk[p] = np.zeros(G)
        if p not in true_pseudobulk:
            warnings.warn(f"Perturbation {p} missing in true data; using zeros vector.")
            true_pseudobulk[p] = np.zeros(G)

    if control_label not in true_pseudobulk or control_label not in pred_pseudobulk:
        raise ValueError("Control label not present in pseudobulk results for true or pred.")

    # Calculate delta for true and predicted perturbations
    true_delta = {p: true_pseudobulk[p] - true_pseudobulk[control_label] for p in perturbations}
    pred_delta = {p: pred_pseudobulk[p] - pred_pseudobulk[control_label] for p in perturbations}

    # Compute cosine similarity for ranking
    def cosine(a, b):
        na = np.linalg.norm(a)
        nb = np.linalg.norm(b)
        if na == 0 or nb == 0:
            return 0
        return np.dot(a, b) / (na * nb)

    PDS_vals = []
    for p in perturbations:
        # Calculate cosine similarity between predicted and true perturbation deltas
        # similarities = [np.dot(pred_delta[p], true_delta[q]) / (np.linalg.norm(pred_delta[p]) * np.linalg.norm(true_delta[q])) for q in perturbations]
        similarities = [cosine(pred_delta[p], true_delta[q]) for q in perturbations]
        # Rank based on similarity
        sim_dict = dict(zip(perturbations, similarities))
        ranked = sorted(perturbations, key=lambda q: sim_dict[q], reverse=True)
        rank = ranked.index(p) + 1

        # Normalize PDS to [0, 1]
        PDS_p = 1.0 - (rank - 1) / float(N)
        PDS_vals.append(PDS_p)

    PDS_array = np.array(PDS_vals, dtype=float)
    PDS_mean = PDS_array.mean()


    # --- 3) 计算 MAE --- (Modifying MAE section)

    # Step 1: 计算所有 perturbations 的 LFC（每个基因取平均绝对 LFC）
    gene_LFC = np.zeros(G)
    for p in perturbations:
        gene_LFC += np.abs(true_delta[p])
    gene_LFC /= len(perturbations)

    # Step 2: 选择 top2000 差异最大的基因
    top_2000_genes = np.argsort(gene_LFC)[-2000:]

    # Step 3: Calculate MAE for each perturbation based on the top 2000 genes
    MAE_vals = []
    for p in perturbations:
        y_true = true_pseudobulk[p]
        y_pred = pred_pseudobulk.get(p, np.zeros(G))
        
        # Filter the top 2000 genes
        y_true_top2000 = y_true[top_2000_genes]
        y_pred_top2000 = y_pred[top_2000_genes]
        
        # Compute MAE
        mae_k = np.mean(np.abs(y_pred_top2000 - y_true_top2000))
        MAE_vals.append(mae_k)

    MAE_array = np.array(MAE_vals, dtype=float)
    MAE_mean = MAE_array.mean()


    # --- baseline handling ---
    DES_baseline = float(baseline_scores.get('DES', 0.0))
    PDS_baseline = float(baseline_scores.get('PDS', 0.0))
    MAE_baseline = float(baseline_scores.get('MAE', 1.0))

    def safe_scale_score(pred, base):
        if base >= 1.0 - 1e-12:
            return 0.0
        return (pred - base) / (1.0 - base)

    DES_scaled = max(0,safe_scale_score(DES_mean, DES_baseline))
    PDS_scaled = max(0,safe_scale_score(PDS_mean, PDS_baseline))
    MAE_scaled = max(0,(MAE_baseline - MAE_mean) / MAE_baseline)

    overall_score = (DES_scaled + PDS_scaled + MAE_scaled) / 3.0 * 100.0

    result = {
        'DES': float(DES_mean),
        'PDS': float(PDS_mean),
        'MAE': float(MAE_mean),
        'DES_scaled': float(DES_scaled),
        'PDS_scaled': float(PDS_scaled),
        'MAE_scaled': float(MAE_scaled),
        'overall_score_percent': float(overall_score),
        'DES_per_perturbation': pd.Series(DES_array, index=perturbations),
        'PDS_per_perturbation': pd.Series(PDS_array, index=perturbations),
        'MAE_per_perturbation': pd.Series(MAE_array, index=perturbations),
        'perturbations': perturbations
    }
    return result

def main_eval(true_file, pred_file, output_path, true_de_csv: Optional[str]=None, pred_de_csv: Optional[str]=None):
    print("Loading AnnData files...")
    baseline_scores = {'DES': 0.0442, 'PDS': 0.5167, 'MAE': 0.1258}

    true_adata = ad.read_h5ad(true_file)
    pred_adata = ad.read_h5ad(pred_file)

    df_true = parallel_differential_expression(
        true_adata,
        groupby_key="target_gene",
        reference="non-targeting",
        threads=8,                 # C++ multithreading
        tie_correction=True,
        use_continuity=True,       # continuity correction for U -> Z
        show_progress=True,
    )
    df_true.to_csv(f"{true_de_csv}", index=False)
    df_pred=parallel_differential_expression(
        pred_adata,
        groupby_key="target_gene",
        reference="non-targeting",
        threads=8,                 # C++ multithreading
        tie_correction=True,
        use_continuity=True,       # continuity correction for U -> Z
        show_progress=True,
    )
    df_pred.to_csv(f"{pred_de_csv}", index=False)

    print("Computing scores (this may take time depending on dataset size)...")
    # 将 true_de_csv / pred_de_csv 传入 compute_vcc_scores（可以为 None）
    scores = compute_vcc_scores(true_adata, pred_adata,
                                baseline_scores=baseline_scores,
                                groupby='target_gene',
                                control_label='non-targeting',
                                fdr=0.05,
                                sample_n=1000,
                                min_cells_per_group=3,
                                true_de_df=true_de_csv,
                                pred_de_df=pred_de_csv)

    # Print summary and save unchanged（与原逻辑一致）
    print("\n===== Summary =====")
    print(f"DES (mean): {scores['DES']:.6f}")
    print(f"PDS (mean): {scores['PDS']:.6f}")
    print(f"MAE (mean): {scores['MAE']:.6f}")
    print(f"DES_scaled: {scores['DES_scaled']:.6f}")
    print(f"PDS_scaled: {scores['PDS_scaled']:.6f}")
    print(f"MAE_scaled: {scores['MAE_scaled']:.6f}")
    print(f"Overall score (percent): {scores['overall_score_percent']:.4f}%")

    out_json = {
        'DES': scores['DES'],
        'PDS': scores['PDS'],
        'MAE': scores['MAE'],
        'DES_scaled': scores['DES_scaled'],
        'PDS_scaled': scores['PDS_scaled'],
        'MAE_scaled': scores['MAE_scaled'],
        'overall_score_percent': scores['overall_score_percent']
    }
    with open(f'{output_path}.json', 'w') as f:
        json.dump(out_json, f, indent=2)

    df_des = scores['DES_per_perturbation'].rename("DES").to_frame()
    df_pds = scores['PDS_per_perturbation'].rename("PDS").to_frame()
    df_mae = scores['MAE_per_perturbation'].rename("MAE").to_frame()
    df_all = pd.concat([df_des, df_pds, df_mae], axis=1)
    df_all.to_csv(f'{output_path}.csv')
    print(f"\nSaved {output_path}.json and {output_path}.csv")

if __name__ == '__main__':
    # main_eval("../vcc_data/eval/small_set_lognorm.h5ad",
    #           "../linear_model/pred/small_pred_lr.h5ad",
    #           "eval_outcome/small_set",
    #           "de_result/small_de_results.csv",
    #           "de_result/small_pred_de_results.csv")
    main_eval("../checkpoints/forth_run20251202/eval_step=step=22000-val_loss=val_loss=0.9084.ckpt/adata_real.h5ad",
              "../checkpoints/forth_run20251202/eval_step=step=22000-val_loss=val_loss=0.9084.ckpt/adata_pred.h5ad",
              "../checkpoints/forth_run20251202/eval_step=step=22000-val_loss=val_loss=0.9084.ckpt/name",
              "../checkpoints/forth_run20251202/eval_step=step=22000-val_loss=val_loss=0.9084.ckpt/de_result/adata_T_de_results.csv",
              "../checkpoints/forth_run20251202/eval_step=step=22000-val_loss=val_loss=0.9084.ckpt/de_result/pred_de_results.csv")
