#!/bin/bash
#SBATCH -J test_1123_pertrubation_statistics9_train_valid_test20251122
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

#在CellPerturb_VCC_Proj6目录下sbatch ./scripts/try_1123_perturbation_statistics9.slurm
cd ./statistics || exit 1

python analyze_perturbations.py -i ../STATE/train_valid_test_split20251122/training_set_1122.h5ad -o train20251122_perturbation_analysis

python analyze_perturbations.py -i ../STATE/train_valid_test_split20251122/validation_set_1122.h5ad -o validation20251122_perturbation_analysis

python analyze_perturbations.py -i ../STATE/train_valid_test_split20251122/test_set_1122.h5ad -o test20251122_perturbation_analysis

# cut -d',' -f1 ./statistics/validation20251122_perturbation_analysis/gene_cell_counts.csv \
  > ./statistics/validation20251122_perturbation_analysis/valid20251122_gene_names.csv

# cut -d',' -f1 ./statistics/test20251122_perturbation_analysis/gene_cell_counts.csv \
  > ./statistics/test20251122_perturbation_analysis/test20251122_gene_names.csv

