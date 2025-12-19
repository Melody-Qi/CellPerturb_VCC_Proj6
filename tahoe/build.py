import os
import glob
import json
import ast
import random
from collections import defaultdict

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from tqdm import tqdm

# ----------------------------
# Config
# ----------------------------
LOCAL_TAHOE_DIR = "./tahoe_small_download"
TRAIN_GLOB = os.path.join(LOCAL_TAHOE_DIR, "data", "train-000[0-9][0-9]-of-03388.parquet")

OUT_DIR = "./tahoe_under_10gb_subset"
os.makedirs(OUT_DIR, exist_ok=True)

TARGET_DISK_GB = 9.5
MAX_PER_CLASS = 500
MAX_CLASSES = 300
TOP_K_GENES = 2048
ROWS_PER_SHARD = 50_000
RANDOM_SEED = 42

KEEP_COLS = [
    "gene_ids",
    "counts",
    "y",
    "target_gene",
    "drug",
    "cell_line_id",
]

parquet_write_kwargs = dict(compression="zstd", compression_level=3)


# ----------------------------
# Helpers
# ----------------------------
def safe_list(x):
    """Normalize targets field to a python list."""
    if x is None:
        return []
    # Tahoe drug_metadata 里 targets 可能直接是单个基因字符串，比如 'PSMB5'
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return []
        # 如果包含常见分隔符，视为多靶/不确定，丢弃
        if any(sep in s for sep in [",", ";", "|", " "]):
            return []
        return [s]
    if isinstance(x, list):
        return x
    # 兼容 numpy/pyarrow 一类带 tolist 的对象
    if hasattr(x, "tolist"):
        try:
            v = x.tolist()
            if isinstance(v, list):
                return v
            return [v] if v is not None else []
        except Exception:
            return []
    return []



def topk_by_counts(gene_ids, counts, k):
    if len(gene_ids) <= k:
        return gene_ids, counts
    idx = sorted(range(len(counts)), key=lambda i: counts[i], reverse=True)[:k]
    return [gene_ids[i] for i in idx], [counts[i] for i in idx]


def get_dir_size_bytes(path: str) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for fn in files:
            fp = os.path.join(root, fn)
            try:
                total += os.path.getsize(fp)
            except OSError:
                pass
    return total


# ----------------------------
# 1) Load local drug_metadata
# ----------------------------
drug_md_path = os.path.join(LOCAL_TAHOE_DIR, "metadata", "drug_metadata.parquet")
if not os.path.exists(drug_md_path):
    raise FileNotFoundError(f"Missing: {drug_md_path}")

print("Loading local drug_metadata ...")
drug_df = pq.read_table(drug_md_path).to_pandas()

drug_df["targets_norm"] = drug_df["targets"].apply(safe_list)

drug_df = drug_df[
    (drug_df["drug"] != "DMSO_TF") &
    (drug_df["targets_norm"].map(len) == 1)
].copy()
drug_df["target_gene"] = drug_df["targets_norm"].map(lambda xs: xs[0])


drug2target = dict(zip(drug_df["drug"], drug_df["target_gene"]))

all_targets = sorted(set(drug2target.values()))
if MAX_CLASSES and len(all_targets) > MAX_CLASSES:
    all_targets = all_targets[:MAX_CLASSES]
target_set = set(all_targets)

drug2target = {d: tg for d, tg in drug2target.items() if tg in target_set}
targets_sorted = sorted(set(drug2target.values()))
gene2y = {g: i for i, g in enumerate(targets_sorted)}
y2gene = {i: g for g, i in gene2y.items()}

print(f"Single-target classes kept: {len(gene2y)}")
if len(gene2y) == 0:
    # Quick debug: show a few raw targets to confirm format
    print("DEBUG: Example raw targets values:")
    print(drug_df["targets"].head(10).tolist())
    raise RuntimeError("No classes found. targets parsing likely still wrong for your local file.")


# ----------------------------
# 2) Scan local parquet shards with pyarrow.dataset
# ----------------------------
train_files = sorted(glob.glob(TRAIN_GLOB))
if not train_files:
    raise FileNotFoundError(f"No files matched TRAIN_GLOB: {TRAIN_GLOB}")

random.seed(RANDOM_SEED)
random.shuffle(train_files)  # cheap shuffle across files

print(f"Found {len(train_files)} local train parquet shards.")
dataset = ds.dataset(train_files, format="parquet")

# only read needed columns
COLUMNS = ["drug", "cell_line_id", "genes", "expressions"]

# ----------------------------
# 3) Sample & write output shards
# ----------------------------
count_per_y = defaultdict(int)
total_kept = 0
buffer = []
shard_idx = 0


def flush():
    global buffer, shard_idx
    if not buffer:
        return
    table = pa.Table.from_pylist(buffer)
    out_path = os.path.join(OUT_DIR, f"train_{shard_idx:04d}.parquet")
    pq.write_table(table, out_path, **parquet_write_kwargs)
    shard_idx += 1
    buffer = []


def all_filled():
    if len(count_per_y) < len(gene2y):
        return False
    return all(count_per_y[i] >= MAX_PER_CLASS for i in range(len(gene2y)))


# Iterate in record batches (efficient)
scanner = dataset.scanner(columns=COLUMNS, batch_size=4096)
print("Sampling cells (class-balanced) from local parquet (pyarrow) ...")

pbar = tqdm(total=None)
for batch in scanner.to_batches():
    # Convert columns to Python lists
    drugs = batch.column(batch.schema.get_field_index("drug")).to_pylist()
    cell_line_ids = batch.column(batch.schema.get_field_index("cell_line_id")).to_pylist()
    genes_col = batch.column(batch.schema.get_field_index("genes")).to_pylist()
    exprs_col = batch.column(batch.schema.get_field_index("expressions")).to_pylist()

    # Shuffle within batch for extra randomness
    order = list(range(len(drugs)))
    random.shuffle(order)

    for i in order:
        drug = drugs[i]
        if drug not in drug2target:
            continue

        target_gene = drug2target[drug]
        y = gene2y[target_gene]
        if count_per_y[y] >= MAX_PER_CLASS:
            continue

        genes = genes_col[i] or []
        exprs = exprs_col[i] or []
        if len(genes) != len(exprs) or len(genes) == 0:
            continue

        # drop CLS/marker
        genes = genes[1:]
        exprs = exprs[1:]
        if len(genes) == 0:
            continue

        genes_k, exprs_k = topk_by_counts(genes, exprs, TOP_K_GENES)

        row = {
            "gene_ids": genes_k,
            "counts": exprs_k,
            "y": int(y),
            "target_gene": target_gene,
            "drug": drug,
            "cell_line_id": cell_line_ids[i],
        }
        row = {k: row[k] for k in KEEP_COLS}

        buffer.append(row)
        count_per_y[y] += 1
        total_kept += 1
        pbar.update(1)

        if len(buffer) >= ROWS_PER_SHARD:
            flush()
            size_gb = get_dir_size_bytes(OUT_DIR) / (1024**3)
            if size_gb >= TARGET_DISK_GB:
                print(f"\nReached disk budget ~{size_gb:.2f} GB >= {TARGET_DISK_GB} GB, stopping early.")
                pbar.close()
                flush()
                print("Done.")
                print(f"Output dir: {OUT_DIR}")
                print(f"Total cells kept: {total_kept:,}")
                print(f"Shards written: {shard_idx}")
                raise SystemExit(0)

        if all_filled():
            pbar.close()
            flush()
            print("\nAll classes filled. Done.")
            print(f"Output dir: {OUT_DIR}")
            print(f"Total cells kept: {total_kept:,}")
            print(f"Shards written: {shard_idx}")
            raise SystemExit(0)

pbar.close()
flush()

print("\nDone (exhausted local files).")
print(f"Output dir: {OUT_DIR}")
print(f"Total cells kept: {total_kept:,}")
print(f"Shards written: {shard_idx}")
print(f"Final disk usage: {get_dir_size_bytes(OUT_DIR)/(1024**3):.2f} GB")

'''
export http_proxy=http://192.168.10.108:10808
export https_proxy=http://192.168.10.108:10808
huggingface-cli login

OUT=./tahoe_small_download
mkdir -p "$OUT"

hf download tahoebio/Tahoe-100M \
  --repo-type dataset \
  --local-dir "$OUT" \
  --include "metadata/drug_metadata.parquet"

hf download tahoebio/Tahoe-100M \
  --repo-type dataset \
  --local-dir "$OUT" \
  --include "data/train-000[1-9][0-9]-of-03388.parquet"
  
hf download tahoebio/Tahoe-100M \
  --repo-type dataset \
  --local-dir "./tahoe_small_download" \
  --include "metadata/gene_metadata.parquet"

result: train-000[1-9][0-9]
Loading local drug_metadata ...
Single-target classes kept: 94
Found 90 local train parquet shards.
Sampling cells (class-balanced) from local parquet (pyarrow) ...
12000it [21:53,  9.14it/s] 

Done (exhausted local files).
Output dir: ./tahoe_under_10gb_subset
Total cells kept: 12,000
Shards written: 1
Final disk usage: 0.03 GB
'''