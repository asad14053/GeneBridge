#!/bin/bash
#SBATCH --job-name=vista_b200_fixed
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

echo "========================================"
echo "Job started on: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-not_set}"
echo "========================================"

ENV_PATH="/scratch/mjabin/conda_envs/Imputation"
PYTHON="${ENV_PATH}/bin/python3.9"

# Clean possible Python contamination
unset PYTHONHOME
unset PYTHONPATH

cd /users/mjabin/projects/GeneBridge

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

echo ""
echo "===== ENVIRONMENT DEBUG ====="
echo "ENV_PATH=${ENV_PATH}"
echo "PYTHON=${PYTHON}"

echo ""
echo "Python files:"
ls -l "${ENV_PATH}/bin/python"*
namei -l "${PYTHON}"

echo ""
echo "Python executable:"
"${PYTHON}" -c "import sys; print(sys.executable)"

echo ""
echo "Python version:"
"${PYTHON}" --version

echo ""
echo "Python sys.path:"
"${PYTHON}" - <<'PY'
import sys
print("\n".join(sys.path[:10]))
PY

echo ""
echo "===== GPU DEBUG ====="
nvidia-smi || true

echo ""
echo "===== TORCH DEBUG ====="
"${PYTHON}" - <<'PY'
import sys
print("python:", sys.executable)

import torch
print("torch:", torch.__version__)
print("torch file:", torch.__file__)
print("torch CUDA build:", torch.version.cuda)
print("cuda available:", torch.cuda.is_available())
print("gpu count:", torch.cuda.device_count())

if torch.cuda.is_available():
    print("gpu name:", torch.cuda.get_device_name(0))
    print("gpu capability:", torch.cuda.get_device_capability(0))
PY

echo ""
echo "===== RUNNING SCRIPT ====="

SCRIPT="$1"
shift

echo "${PYTHON} ${SCRIPT} $*"

"${PYTHON}" "${SCRIPT}" "$@"

echo ""
echo "Python script finished."
echo "Date: $(date)"