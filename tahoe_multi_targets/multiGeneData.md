# 数据处理流程说明



## 1. 输入数据与关键文件

### 1.1 Tahoe 训练数据（Parquet）
- 位置：`./tahoe_small_download/data/train*.parquet`
- 读取列（`COLUMNS`）：
  - `genes`：基因 token id 列表（List[int]）
  - `expressions`：表达值列表（List[float]）
  - `drug`：药物名（string）
  - `cell_line_id`：细胞系 ID
  - `sample`：样本信息
  - `plate`：板号/批次


### 1.2 药物元数据（drug_metadata.parquet）
- 位置：`./tahoe_small_download/metadata/drug_metadata.parquet`
- 用途：
  - 从 `drug` 映射到 `targets`（靶基因集合）
  - 支持单靶点和多靶点药物（本脚本保留“targets 数量 >= 1”的所有药物）

### 1.3 Tahoe token id → scGPT vocab id 映射表（JSON）
- 文件：`tahoe_tokenid_to_scgptid.json`
- 用途：把输入的 Tahoe `genes` token id 映射成 scGPT 体系的 gene id
- 未出现在映射表中的 gene token 会被丢弃（drop unknown）

---

## 2. 输出目录与输出文件

### 2.1 输出目录
- `OUT_DIR = "./tahoe_mulGene_scgptid_firstDMSO_ratio02"`
- 脚本会自动创建该目录。

### 2.2 主要输出：分片 Parquet 数据集
- 文件名模式：`single_target_0000.parquet`, `single_target_0001.parquet`, ...
- 分片策略：
  - 每个 shard 最多写入 `ROWS_PER_SHARD = 100000` 行
  - 达到阈值触发 `flush()` 写盘

### 2.3 辅助输出：标签词表与 ID 空间声明
- `label_vocab.json`
  - `gene2label`: `{target_gene(str): label(int)}`
  - `label2gene`: `{label(int): target_gene(str)}`
  - 标签空间是所有药物 targets 的并集（union）
- `id_space.json`
  - `{"genes": "scgpt"}` 表示 `genes` 字段已经处于 scGPT ID 空间

---


## 3

- 仅保留 `targets_list` 长度 >= 1 的药物
  - `targets_list` 中每个 target 也会强制转为字符串并移除空字符串
  - 结果：`drug2targets[drug] = [t1, t2, ...]`
- **可选：丢弃第一个 token**
   - `DROP_FIRST_TOKEN = False`（默认不丢）
- **按 gene id 升序排序**
   - 排序基因列表，并同步重排表达值列表
- 规则：同一个 key 只记录第一次遇到的 DMSO 行，后续 DMSO 行直接跳过。
- `overlap_ratio = |tg ∩ cg| / |tg|`
  - 若 `overlap_ratio < OVERLAP_RATIO_MIN (0.2)` → reject
- control 表达对齐到 treat 基因空间
  - `ctrl_genes` 会被直接设为 `treat_genes`（对齐后的基因空间）
  - 若 control 中缺失某个 treat gene，则用 `CTRL_FILL_VALUE = 0.0` 填充
- 多靶点扩增为多行: 对于某条 treat，如果 `targets = [g1, g2, g3]`：
  - 会写出 3 行记录
  - 每行除了 `target_gene` 和 `label` 不同，其余字段相同
- 每个 target_gene 的行数上限 `MAX_ROWS_PER_TARGET = 100000`
  - 若某个 `target_gene` 已累计写出行数达到上限，则后续同靶点行会被拒绝（`cap_rejected += 1`）

---

## 最终输出数据格式（Parquet schema 级别）

每行输出代表一个 **(treat, control) paired 样本**，并且已经展开到“单 target_gene”粒度。

输出列（`KEEP_COLS`）如下：

###  treat 相关字段
- `genes`: `List[int]`  
  - treat 的 scGPT gene ids（升序）
- `expressions`: `List[float]`  
  - treat 表达值（与 genes 对齐）
- `drug`: `string`  
  - treat 药物名称
- `target_gene`: `string`  
  - 当前这一行对应的靶基因（多靶点药物会多行）
- `label`: `int`  
  - `target_gene` 对应的类别 id（来自 `label_vocab.json`）
- `cell_line_id`: `string | null`
- `sample`: string
- `plate`: string

###  control 相关字段（已对齐到 treat 基因空间）
- `ctrl_genes`: `List[int]`
  - 与 `genes` 相同（对齐后的基因空间）
- `ctrl_expressions`: `List[float]`
  - control 表达值，按 `genes` 顺序对齐，不存在则填 0.0
- `ctrl_drug`: `string`
  - control 药物名称（通常含 DMSO）

### 过滤指标字段
- `overlap_ratio`: `float`
  - $ |treat \cap ctrl| / |treat| $

---

## 你可以如何使用输出数据

- 训练分类模型：
  - 输入：`genes`, `expressions`, `ctrl_expressions`（或 treat/control 差分特征）
  - 标签：`label`
- 也可以按 `target_gene` 分组统计（注意已做了 target 级展开）
- 多靶点药物会重复样本（同一个 treat-control pair）到多个 target 上，这是设计使然

---

## 关键配置项速览

- `OVERLAP_RATIO_MIN = 0.2`：treat-control 基因重叠比过滤阈值
- `ROWS_PER_SHARD = 100000`：输出分片行数
- `MAX_ROWS_PER_TARGET = 100000`：每个 target_gene 的最大保留行数
- 默认：`key = (cell_line_id, plate)`  
- `CTRL_FILL_VALUE = 0.0`：对齐时 control 缺失基因的填充值
- `DROP_FIRST_TOKEN = False`：是否丢弃 gene 列表首 token（默认否）

---
