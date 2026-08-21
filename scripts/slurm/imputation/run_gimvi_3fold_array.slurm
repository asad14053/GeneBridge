#!/bin/bash
#SBATCH --job-name=gimvi_3fold
#SBATCH --output=logs/gimvi_%A_%a.out
#SBATCH --error=logs/gimvi_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --array=1-3%3

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

RUNNER="${PROJECT_ROOT}/src/imputation/gimVI/run_gimvi_fold.py"

EXPERIMENT="${1:?Missing experiment key}"
EXPERIMENT_LABEL="${2:?Missing experiment label}"
REFERENCE="${3:?Missing reference H5AD}"
RUN_MODE="${4:-full}"

FOLD="${SLURM_ARRAY_TASK_ID}"

export PYTHONNOUSERSITE=1
export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

cd "${PROJECT_ROOT}"

mkdir -p logs

echo "===================================================================================================="
echo "gimVI three-fold benchmark"
echo "===================================================================================================="
echo "Array job: ${SLURM_ARRAY_JOB_ID}"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo
echo "Experiment: ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"
echo "Fold: ${FOLD}"
echo "Reference: ${REFERENCE}"
echo "Run mode: ${RUN_MODE}"
echo "Memory: 200G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Array policy: 1-3%3"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${RUNNER}" || {
    echo "ERROR: gimVI runner not found:"
    echo "${RUNNER}"
    exit 1
}

test -f "${REFERENCE}" || {
    echo "ERROR: Reference H5AD not found:"
    echo "${REFERENCE}"
    exit 1
}

"${PYTHON}" - <<'PY'
import scvi
import torch
from scvi.external import GIMVI

print("scvi-tools:", scvi.__version__)
print("PyTorch:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())
print("GIMVI import: PASS")
PY

"${PYTHON}" -m py_compile \
  "${RUNNER}"

EXTRA_ARGS=()

if [[ "${RUN_MODE}" == "smoke" ]]
then
    EXTRA_ARGS+=(--skip-model-save)
fi

srun --cpu-bind=cores \
  "${PYTHON}" -u "${RUNNER}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --experiment-label "${EXPERIMENT_LABEL}" \
    --fold "${FOLD}" \
    --reference "${REFERENCE}" \
    --reference-mode benchmark_300 \
    --run-mode "${RUN_MODE}" \
    --smoke-epochs 2 \
    --full-epochs 200 \
    --batch-size 128 \
    --n-latent 32 \
    --train-size 0.90 \
    --learning-rate 1e-4 \
    --gradient-clip-val 5.0 \
    --seed 8667 \
    --cpus 16 \
    --k-spatial 15 \
    --ssim-grid-size 128 \
    --n-clusters 10 \
    --n-pcs 30 \
    --plot-genes 10 \
    "${EXTRA_ARGS[@]}"

echo
echo "Completed: $(date)"
