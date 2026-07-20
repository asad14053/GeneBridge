#!/bin/bash
#SBATCH --job-name=ipynb_cpu
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail
mkdir -p logs outputs/notebooks

source ~/anaconda3/etc/profile.d/conda.sh
conda activate /scratch/mjabin/conda_envs/Imputation

cd /users/mjabin/projects/GeneBridge
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

NOTEBOOK="$1"
BASENAME=$(basename "$NOTEBOOK" .ipynb)

jupyter nbconvert \
  --to notebook \
  --execute "$NOTEBOOK" \
  --output "${BASENAME}_EXECUTED.ipynb" \
  --output-dir outputs/notebooks \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.kernel_name=Imputation
