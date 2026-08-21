#!/bin/bash
#SBATCH --job-name=envi_3fold
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-3%3

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage:"
    echo "sbatch --job-name=envi_ex5 $0 EXPERIMENT EXPERIMENT_LABEL REFERENCE"
    exit 1
fi

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

SCRIPT="${PROJECT_ROOT}/src/imputation/ENVI/run_envi_fold.py"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"

OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

EXPERIMENT="$1"
EXPERIMENT_LABEL="$2"
REFERENCE="$3"

FOLD="${SLURM_ARRAY_TASK_ID}"

cd "${PROJECT_ROOT}"

mkdir -p \
  logs \
  "${OUTPUT_ROOT}"


# ============================================================
# CPU / JAX
# ============================================================

export JAX_PLATFORMS=cpu
export JAX_PLATFORM_NAME=cpu
export XLA_PYTHON_CLIENT_PREALLOCATE=false

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1


# ============================================================
# Automatically find working ENVI environment
# ============================================================

PYTHON=""

CANDIDATES=(
    "/beegfs/labs/hulab/projects/mjabin/conda_envs/envi/bin/python"
    "/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"
)

for CANDIDATE in "${CANDIDATES[@]}"; do

    if [[ -x "${CANDIDATE}" ]] && \
       "${CANDIDATE}" -c \
       'import scenvi, jax, scanpy, anndata; assert jax.default_backend().lower() == "cpu"' \
       >/dev/null 2>&1
    then
        PYTHON="${CANDIDATE}"
        break
    fi

done


if [[ -z "${PYTHON}" ]]; then

    echo "ERROR: No working ENVI Python environment found."
    echo
    echo "Checked:"

    printf '  %s\n' "${CANDIDATES[@]}"

    echo
    echo "Required imports:"
    echo "  scenvi"
    echo "  jax"
    echo "  scanpy"
    echo "  anndata"

    exit 1
fi


# ============================================================
# Log configuration
# ============================================================

echo "===================================================================================================="
echo "ENVI three-fold benchmark"
echo "===================================================================================================="

echo "Job: ${SLURM_JOB_ID}"
echo "Array job: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Task/fold: ${FOLD}"
echo "Node: $(hostname)"
echo "Started: $(date)"

echo
echo "Experiment: ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"
echo "Reference: ${REFERENCE}"
echo "Python: ${PYTHON}"

echo
echo "Training steps: 16000"
echo "Batch size: 128"
echo "num_HVG: 2048"
echo "num_cov_genes: 64"
echo "k_nearest: 8"
echo "covet_batch_size: 256"
echo "spatial_dist: pois"
echo "sc_dist: nb"
echo "stable_eps: 1e-6"

echo
echo "JAX backend: cpu"
echo "Memory: 200G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Array policy: 1-3%3"

echo "===================================================================================================="


# ============================================================
# Run
# ============================================================

"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --experiment-label "${EXPERIMENT_LABEL}" \
    --fold "${FOLD}" \
    --reference "${REFERENCE}" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --run-mode full \
    --full-training-steps 16000 \
    --batch-size 128 \
    --num-hvg 2048 \
    --num-cov-genes 64 \
    --k-nearest 8 \
    --covet-batch-size 256 \
    --spatial-dist pois \
    --sc-dist nb \
    --stable-eps 1e-6 \
    --seed 8667 \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --k-spatial 15 \
    --ssim-grid-size 128 \
    --n-clusters 10 \
    --n-pcs 30 \
    --plot-genes 10


echo
echo "===================================================================================================="
echo "SUCCESS"
echo "Experiment: ${EXPERIMENT}"
echo "Fold: ${FOLD}"
echo "Completed: $(date)"
echo "===================================================================================================="
