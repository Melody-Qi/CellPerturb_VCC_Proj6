import os
import glob
import pyarrow.parquet as pq

# 要处理的目录列表
INPUT_DIRS = [
    ".",
    "../tahoe_single_target_scgptid_firstDMSO_ratio05",
]

def convert_parquet_to_csv(parquet_path):
    parquet_path = os.path.abspath(parquet_path)
    parent_dir = os.path.dirname(parquet_path)

    # 输出到 parquet 所在目录下的 csv/ 子目录
    csv_dir = os.path.join(parent_dir, "csv")
    os.makedirs(csv_dir, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(parquet_path))[0]
    csv_path = os.path.join(csv_dir, base_name + ".csv")

    print(f"Converting: {parquet_path}")
    print(f" -> {csv_path}")

    # 读取 parquet 并写 csv
    table = pq.read_table(parquet_path)
    table.to_pandas().to_csv(csv_path, index=False)

def main():
    for input_dir in INPUT_DIRS:
        input_dir = os.path.abspath(input_dir)
        print(f"\nScanning directory: {input_dir}")

        parquet_files = glob.glob(
            os.path.join(input_dir, "**", "*.parquet"),
            recursive=True
        )

        print(f"Found {len(parquet_files)} parquet files")

        for parquet_file in parquet_files:
            convert_parquet_to_csv(parquet_file)

if __name__ == "__main__":
    main()
