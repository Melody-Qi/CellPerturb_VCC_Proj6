import glob
import gc
import pickle
from collections import defaultdict
import pandas as pd
import os

PARQUET_FILES = sorted(glob.glob("single_target_*.parquet"))
GENE_META_CSV = "gene_metadata.csv"
CACHE_DIR = "parquet_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

# ===============================
# load gene metadata
# ===============================
gene_meta = pd.read_csv(GENE_META_CSV)
symbol2token = dict(zip(gene_meta["gene_symbol"], gene_meta["token_id"]))
print(f"Loaded {len(symbol2token)} gene symbols")

# ===============================
# process each parquet independently
# ===============================
for f in PARQUET_FILES:
    print(f"\nProcessing {f}")
    df = pd.read_parquet(
        f,
        columns=["genes", "expressions", "ctrl_genes", "ctrl_expressions", "drug", "target_gene", "cell_line_id"]
    )

    # per-parquet stats
    stats = {}

    for _, row in df.iterrows():
        cell, drug = row["cell_line_id"], row["drug"]
        key = (cell, drug)

        target_symbol = row["target_gene"]
        target_token = symbol2token.get(target_symbol)
        if target_token is None:
            print("What? ", drug,"的 target gene 没找到？")
            continue

        if key not in stats:
            stats[key] = {
                "pert_sum": defaultdict(float),
                "pert_cnt": defaultdict(int),
                "ctrl_sum": defaultdict(float),
                "ctrl_cnt": defaultdict(int),
                "target_token": target_token,
                "target_symbol": target_symbol,
                "row_count": 0,
            }

        s = stats[key]
        s["row_count"] += 1
        for g, x in zip(row["genes"], row["expressions"]):
            g = int(g)
            s["pert_sum"][g] += x
            s["pert_cnt"][g] += 1
        for g, x in zip(row["ctrl_genes"], row["ctrl_expressions"]):
            g = int(g)
            s["ctrl_sum"][g] += x
            s["ctrl_cnt"][g] += 1

    del df
    gc.collect()

    # compute means/deltas for this parquet
    parquet_means = {}
    for key, s in stats.items():
        pert_mean = {g: s["pert_sum"][g] / s["pert_cnt"][g] for g in s["pert_sum"]}
        ctrl_mean = {g: s["ctrl_sum"][g] / s["ctrl_cnt"][g] for g in s["ctrl_sum"]}
        all_genes = set(pert_mean) | set(ctrl_mean)
        delta = {g: pert_mean.get(g, 0.0) - ctrl_mean.get(g, 0.0) for g in all_genes}

        parquet_means[key] = {
            "pert_mean": pert_mean,
            "ctrl_mean": ctrl_mean,
            "delta": delta,
            "target_token": s["target_token"],
            "target_symbol": s["target_symbol"],
            "row_count": s["row_count"],
        }

    # save this parquet's summary
    parquet_name = os.path.basename(f).replace(".parquet", ".pkl")
    cache_file = os.path.join(CACHE_DIR, parquet_name)
    with open(cache_file, "wb") as fout:
        pickle.dump(parquet_means, fout)

    print(f"✅ Saved {len(parquet_means)} keys for {f} to {cache_file}")

    # release memory
    del stats, parquet_means
    gc.collect()

print("\nAll parquet files processed individually. Each saved in cache directory.")


