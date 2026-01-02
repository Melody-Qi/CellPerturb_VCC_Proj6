#!/usr/bin/env python
"""
Aggregate all parquet_cache summaries into global (cell_line_id, drug) statistics.

Weighted by each parquet's row_count for accurate averaging.
"""

import os
import glob
import pickle
from collections import defaultdict

CACHE_DIR = "parquet_cache"

# -----------------------------
# load all cached parquet summaries
# -----------------------------
cache_files = sorted(glob.glob(os.path.join(CACHE_DIR, "*.pkl")))
print(f"Found {len(cache_files)} cached parquet summaries")

# structure: agg_stats[(cell_line_id, drug)] = {"pert_sum", "ctrl_sum", "total_count", "target_token", "target_symbol"}
agg_stats = {}

for cache_file in cache_files:
    with open(cache_file, "rb") as f:
        parquet_means = pickle.load(f)
    
    for key, s in parquet_means.items():  # key = (cell_line_id, drug)
        cell_line, drug = key
        row_count = s["row_count"]
        
        if key not in agg_stats:
            agg_stats[key] = {
                "pert_sum": defaultdict(float),
                "ctrl_sum": defaultdict(float),
                "total_count": 0
            }
        
        agg = agg_stats[key]
        agg["total_count"] += row_count
        
        # weighted sum
        for g, val in s["pert_mean"].items():
            agg["pert_sum"][g] += val * row_count
        for g, val in s["ctrl_mean"].items():
            agg["ctrl_sum"][g] += val * row_count

print("✅ All cached parquet files loaded and summed")

# -----------------------------
# compute weighted average
# -----------------------------
global_stats = {}

for key, agg in agg_stats.items():
    total_count = agg["total_count"]
    
    # weighted average
    pert_mean = {g: val / total_count for g, val in agg["pert_sum"].items()}
    ctrl_mean = {g: val / total_count for g, val in agg["ctrl_sum"].items()}
    
    all_genes = set(pert_mean) | set(ctrl_mean)
    delta = {g: pert_mean.get(g, 0.0) - ctrl_mean.get(g, 0.0) for g in all_genes}
    
    global_stats[key] = {
        "pert_mean": pert_mean,
        "ctrl_mean": ctrl_mean,
        "delta": delta,
        "total_count": total_count
    }

print(f"✅ Computed global statistics for {len(global_stats)} (cell_line, drug) combinations")

# -----------------------------
# save aggregated results
# -----------------------------
output_file = os.path.join(CACHE_DIR, "global_parquet_stats.pkl")
with open(output_file, "wb") as f:
    pickle.dump(global_stats, f)

print(f"✅ Saved global stats to {output_file}")
