#!/bin/bash
#SBATCH --job-name=vista_gpu
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=rp6b-1-gm96-c8-m64
#SBATCH --account=ai-gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G

set -euo pipefail

mkdir -p logs
mkdir -p outputs/notebooks

echo "Job started on: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-not_set}"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate /scratch/mjabin/conda_envs/Imputation

cd /users/mjabin/projects/GeneBridge

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

NOTEBOOK="$1"
BASENAME=$(basename "$NOTEBOOK" .ipynb)

echo "Running notebook: $NOTEBOOK"

python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
print("gpu count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("gpu name:", torch.cuda.get_device_name(0))
PY

jupyter nbconvert \
  --to notebook \
  --execute "$NOTEBOOK" \
  --output "${BASENAME}_EXECUTED.ipynb" \
  --output-dir outputs/notebooks \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.kernel_name=Imputation

echo "Notebook finished:"
echo "outputs/notebooks/${BASENAME}_EXECUTED.ipynb"
echo "Date: $(date)"
