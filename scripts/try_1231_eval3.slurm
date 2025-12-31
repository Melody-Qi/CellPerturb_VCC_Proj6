#!/bin/bash
#SBATCH -J test_1231_eval3
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

#在CellPerturb_VCC_Proj6目录下sbatch ./scripts/try_1231_eval3.slurm
cd ../sea000123/gene-perturbation-prediction || exit 1

export PYTHONHASHSEED=0

python src/main_new.py \
  --config src/configs/h5ad.yaml \
  --eval_mode h5ad \
  --truth_h5ad /public/home/CS286/qiyt2023-CS286/CellPerturb_VCC_Proj6_new4/STATE/checkpoints/forth_run20251202/eval_step=step=22000-val_loss=val_loss=0.9084.ckpt/adata_real.h5ad \
  --pred_h5ad /public/home/CS286/qiyt2023-CS286/CellPerturb_VCC_Proj6_new4/STATE/checkpoints/forth_run20251202/eval_step=step=22000-val_loss=val_loss=0.9084.ckpt/adata_pred.h5ad \
  --pert_col target_gene \
  --control_label non-targeting

