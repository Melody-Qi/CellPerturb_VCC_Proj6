#!/usr/bin/env python
"""
Fix invalid var schema in existing Tahoe h5ad shards and merge them on disk.
"""

from pathlib import Path
import argparse
import scanpy as sc
import pandas as pd
from anndata.experimental import concat_on_disk


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        default="tahoe_mulGene_scgptid_firstDMSO_ratio02_52drug/h5ad",
        help="Directory with original h5ad shards",
    )
    parser.add_argument(
        "--output",
        default="tahoe_mulGene_scgptid_firstDMSO_ratio02_52drug/h5ad/tahoe.h5ad",
        help="Merged output h5ad",
    )
    parser.add_argument(
        "--max-loaded-elems",
        type=int,
        default=20_000_000,
    )
    return parser.parse_args()


def fix_one_h5ad(in_path: Path, out_path: Path):
    adata = sc.read_h5ad(in_path)

    # --- 修复 var ---
    if "gene_idx" in adata.var.columns:
        gene_ids = adata.var["gene_idx"].astype(str).tolist()
        adata.var = pd.DataFrame(index=gene_ids)
    adata.var.index = adata.var.index.astype(str)

    # --- 修复 sparse matrix indptr ---
    import scipy.sparse as sp

    if sp.issparse(adata.X):
        adata.X = adata.X.asformat('csr')
        adata.X.indptr = adata.X.indptr.astype('int64')

    for layer_name in adata.layers.keys():
        if sp.issparse(adata.layers[layer_name]):
            adata.layers[layer_name] = adata.layers[layer_name].asformat('csr')
            adata.layers[layer_name].indptr = adata.layers[layer_name].indptr.astype('int64')

    adata.write_h5ad(out_path)
    del adata



def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fixed_dir = input_dir / "fixed"
    fixed_dir.mkdir(exist_ok=True)

    fixed_paths = []

    print("=== Fixing individual h5ad files ===")
    for p in sorted(input_dir.glob("single_target_*.h5ad")):
        out_p = fixed_dir / p.name
        print(f"Fixing {p.name}")
        fix_one_h5ad(p, out_p)
        fixed_paths.append(str(out_p))

    print("\n=== Merging fixed h5ad files ===")
    concat_on_disk(
        fixed_paths,
        str(output_path),
        axis=0,
        join="inner",
        merge="same",
        uns_merge="same",
        max_loaded_elems=args.max_loaded_elems,
    )

    print(f"\n✅ Done. Merged file saved to:\n{output_path}")


if __name__ == "__main__":
    main()
