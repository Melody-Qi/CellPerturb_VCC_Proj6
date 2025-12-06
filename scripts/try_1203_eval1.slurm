#!/bin/bash
#SBATCH -J test_1203_eval1
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

#在CellPerturb_VCC_Proj6目录下sbatch ./scripts/try_1203_eval1.slurm
cd ./STATE/train_valid_test_split20251201 || exit 1

python my_cell_eval3.py

