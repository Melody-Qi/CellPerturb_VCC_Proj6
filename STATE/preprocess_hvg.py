import argparse
import anndata as ad
import scanpy as sc

def main():
    parser = argparse.ArgumentParser(description="Select HVGs and save to new AnnData file.")
    parser.add_argument("--input", type=str, required=True, help="Path to input .h5ad")
    parser.add_argument("--output", type=str, required=True, help="Path to output .h5ad")
    parser.add_argument("--num_hvgs", type=int, default=2000, help="Number of HVGs to select")
    args = parser.parse_args()

    print(f"📥 Loading input AnnData from: {args.input}")
    adata = ad.read_h5ad(args.input)

    print(f"🔍 Selecting top {args.num_hvgs} highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.num_hvgs)

    print("📦 Extracting HVG expression matrix to obsm['X_hvg'] ...")
    adata.obsm["X_hvg"] = adata[:, adata.var.highly_variable].X.toarray()

    print(f"💾 Saving output to: {args.output}")
    adata.write_h5ad(args.output)

    print("✅ Done! HVG preprocessing completed.")

if __name__ == "__main__":
    main()
