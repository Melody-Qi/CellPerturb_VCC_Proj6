#!/bin/bash
#SBATCH -J test_1111_train1
#SBATCH -p CS286
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=all
#SBATCH --mail-user=2943895132@qq.com
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH --gres=gpu:1
#SBATCH --output=./scripts/%x.out
#SBATCH --error=./scripts/%x.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CellPerturb_VCC_Proj6_STATE2

nvidia-smi

echo "=== SLURM环境变量 ==="
echo "SLURM_GPUS_PER_NODE: $SLURM_GPUS_PER_NODE"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

echo "=== PyTorch检测 ==="
python -c "
import torch
print(f'PyTorch检测GPU数: {torch.cuda.device_count()}')
print(f'可用GPU: {[torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]}')
"
#在CellPerturb_VCC_Proj6目录下sbatch ./scripts/try_1111_train1.slurm
cd ./STATE || exit 1

mkdir -p ./checkpoints

state tx train \
  data.kwargs.toml_config_path="examples/fewshot.toml" \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.num_workers=8 \
  data.kwargs.batch_col=batch \
  data.kwargs.pert_col=target_gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.control_pert=non-targeting \
  training.max_steps=20000 \
  training.ckpt_every_n_steps=100 \
  training.batch_size=8 \
  training.lr=1e-4 \
  model.kwargs.hidden_dim=256 \
  model.kwargs.cell_set_len=64 \
  model=pertsets \
  wandb.tags="[vcc_fewshot]" \
  output_dir="./checkpoints" \
  name="vcc_fewshot"

# state tx train \
# data.kwargs.toml_config_path="examples/fewshot.toml" \
# data.kwargs.embed_key=X_hvg \   指定在 .h5ad 文件中要用的表达矩阵 key，一般是 高变基因矩阵（highly variable genes），例如 adata.obsm["X_hvg"]
# data.kwargs.num_workers=12 \
# data.kwargs.batch_col=batch_var \指定 .h5ad.obs 中的列名，用来区分 实验批次。有时可以为空（如果你没有 batch 信息就不要加）。
# data.kwargs.pert_col=target_gene哪一列表示扰动基因 \
# data.kwargs.cell_type_key=cell_type哪一列表示细胞类型 \
# data.kwargs.control_pert=TARGET1 哪个扰动代表“对照组”