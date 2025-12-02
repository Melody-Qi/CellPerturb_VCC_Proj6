#!/bin/bash
#SBATCH -J test_1202_train5
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
#在CellPerturb_VCC_Proj6目录下sbatch ./scripts/try_1202_train5.slurm
cd ./STATE || exit 1

state tx train \
  data.kwargs.toml_config_path="train_valid_test_split20251201/fewshot.toml" \
  data.kwargs.num_workers=12 \
  data.kwargs.batch_col="batch_var" \
  data.kwargs.pert_col="target_gene" \
  data.kwargs.cell_type_key="cell_type" \
  data.kwargs.control_pert="non-targeting" \
  data.kwargs.perturbation_features_file="competition_support_set/ESM2_pert_features.pt" \
  training.max_steps=60000 \
  training.val_freq=100 \
  training.ckpt_every_n_steps=100 \
  training.lr=1e-4 \
  model=state_sm \
  wandb.tags="[forth_run20251202]" \
  wandb.project=vcc \
  wandb.entity=arcinstitute \
  output_dir="./checkpoints" \
  name="forth_run20251202"

# state tx train \
#   model_dir="./ST-Tahoe" \
#   data.kwargs.toml_config_path="examples/fewshot.toml" \
#   data.kwargs.embed_key=X_hvg \
#   data.kwargs.batch_col=batch \
#   data.kwargs.pert_col=target_gene \
#   data.kwargs.cell_type_key=cell_type \
#   data.kwargs.control_pert=CONTROL \
#   training.max_steps=5000 \
#   training.batch_size=4 \
#   output_dir="./finetune_output" \
#   name="finetune_tahoe"