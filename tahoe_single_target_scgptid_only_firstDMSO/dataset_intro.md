下面按“数据从哪来 → 怎么被过滤/映射 → 怎么配对 → 怎么写出”的顺序，把这段脚本的**完整数据处理流程**讲清楚，并指出每个关键步骤在做什么、为什么这样做、以及会丢掉哪些数据。

---

## 1) 输入与整体目标

### 输入数据

* 本地 Tahoe 数据目录：`./tahoe_small_download`
* 训练数据 parquet：`./tahoe_small_download/data/train*.parquet`
* 药物元数据 parquet：`./tahoe_small_download/metadata/drug_metadata.parquet`
* 基因 ID 映射表 JSON：`tahoe_tokenid_to_scgptid.json`
  * 表示 **Tahoe 的 token_id → scGPT vocab gene id** 的映射

### 输出目标

输出一个新的 parquet 数据集到：
`./tahoe_single_target_scgptid_only_firstDMSO`

每个样本包含：
* treat 的基因列表与表达
* 与 treat 基因顺序对齐后的 ctrl 表达（缺失补 0）
* drug / target_gene / label 等监督信息
* overlap_ratio 作为质量过滤指标
  * ctrl组和药物扰动基因的overlap_ratio小于0.05的数据将会被剔除



## 2) 预处理：构造“单靶点药物”监督空间

- 只保留“单靶点药物”
- 找不到ctrl组的数据会被丢弃
- 给 target_gene 编码成分类标签
- 基因本来对应的是tahoe的id，改成与scgpt的id对应
  * `label_vocab.json`
  * `id_space.json`（声明 genes 使用 scgpt id 空间）
  - 如果有基因在scgpt的id space中未出现，则丢弃该基因（因此会改变原始基因集合大小）
- 可选 `DROP_FIRST_TOKEN`：把第一个 token 丢掉（默认 False）
- 按 scgpt gene id 排序：`order = np.argsort(mapped)`

之前的想法是每个对照组采样200个数据，从里面把扰动后的每一个基因的对照组表达量都提出来，组成一个对照；但是这个太慢了
后来改了一下，直接使用读取到的第一个对照，对照里的基因列表和扰动基因列表重合度挺低的；我现在采用的算法overlap ratio=0.05，没有数据被丢弃；如果设置成0.2，会有一半以上数据不符合要求。

## 3) 为每个 key 收集第一条 DMSO 对照

- `is_control_drug(drug)`：只要 drug 字符串 upper 后包含 `"DMSO"` 就算对照。
-`make_key(cell_line_id, plate, sample)`：
  * 默认 `USE_SAMPLE_IN_KEY = False` → key 为 `(cell_line_id, plate)`
  * 这个 key 决定“treat 要配哪个 ctrl”。
- 为了加快处理速度，这里选择在数据中扫描到的**第一条**key相符的DMSO作为对照
  - 因为储存所有key相符的DMSO，然后统一处理特别慢
  - 会丢弃不合法的数据
  * `ctrl_first`: dict
  * `key -> (ctrl_genes_sorted_scgpt, ctrl_exprs_sorted, ctrl_drug_str, hash)`


## 4) 
### 4.1 overlap_ratio 质量过滤（避免 ctrl/treat 基因空间差太大）
treat 的基因里至少要有 5% 也出现在 ctrl 里，否则 ctrl 对齐会大量补 0，信息太差。

### 4.2 ctrl 表达对齐到 treat 基因顺序

`align_ctrl_to_treat_sorted(treat_genes=tg, ctrl_genes=cg, ctrl_exprs=cx, fill_value=0.0)`：

双指针扫描：

* 如果 treat gene == ctrl gene → 放入 ctrl_expr
* 如果 treat gene < ctrl gene → treat 里多出来的 gene，ctrl 没有 → 填 0
* 如果 treat gene > ctrl gene → ctrl gene 在 treat 不存在 → ctrl 指针前进（等同丢弃）

输出：
* `aligned_genes` 直接等于 treat 的 `tg`
* `aligned_ctrl_exprs` 长度与 tg 相同，对齐后的 ctrl 表达

> 这一步让模型输入可以是：
> treat: (tg, tx)
> ctrl:  (tg, aligned_ctrl_exprs)
> 两者基因轴完全一致。

### 4.3 生成监督 label 并写入缓冲

* `label = gene2label[target_gene]`

写入列缓冲 `colbuf`（columnar），包含：

* treat: genes/expressions/drug/target_gene/label/cell_line_id/sample/plate
* ctrl: ctrl_genes(=tg)/ctrl_expressions(aligned)/ctrl_drug
* overlap_ratio

### 4.4 分片写出 parquet

当缓冲行数达到 `ROWS_PER_SHARD=500_000` 就 `flush()`：
* `pa.Table.from_pydict(colbuf)` → Arrow Table
* `pq.write_table(..., compression="zstd", compression_level=3)`
* 输出文件名：`single_target_0000.parquet`, `single_target_0001.parquet`, ...

最后 loop 结束再 flush 一次，打印统计：
* `kept`=2549364
* `rejected`（只统计 overlap_ratio filter 拒绝的数量）=0

