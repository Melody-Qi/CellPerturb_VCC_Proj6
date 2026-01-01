#!/bin/bash
#SBATCH -J test_0101_topKhit1
#SBATCH -p CS286
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
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

#在CellPerturb_VCC_Proj6目录下sbatch ./scripts/try_0101_topKhit1.slurm
cd ./tahoe_multi_targets/tahoe_mulGene_scgptid_firstDMSO_ratio02_52drug || exit 1

python topK_hit.py

