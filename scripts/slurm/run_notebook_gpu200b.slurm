#!/bin/bash
#SBATCH --job-name=vista_py_gpu
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=b200-8-gm1432-c192-m2048
#SBATCH --account=ai-gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

mkdir -p logs

echo "Job started on: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-not_set}"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate /scratch/mjabin/conda_envs/Imputation

cd /users/mjabin/projects/GeneBridge

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

SCRIPT="$1"
shift

echo "Running Python script:"
echo "python ${SCRIPT} $@"

python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
print("gpu count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("gpu name:", torch.cuda.get_device_name(0))
PY

python "$SCRIPT" "$@"

echo "Python script finished."
echo "Date: $(date)"