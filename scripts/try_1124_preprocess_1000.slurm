#!/bin/bash
#SBATCH -J test_1124_preprocess_1000
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
#在CellPerturb_VCC_Proj6目录下sbatch ./scripts/try_1124_preprocess_1000.slurm
cd ./STATE || exit 1
state tx preprocess_train \
  --adata ./vcc_data/adata_Training.h5ad \
  --output ./preprocessed_data/preprocessed_training_data_1000.h5ad \
  --num_hvgs 1000

# 这一步会做：

# Normalize (sc.pp.normalize_total)
# Log1p (Log transform)
# 选取高变基因 (sc.pp.highly_variable_genes) 1000
# 输出 .obsm['X_hvg']，供模型训练使用。