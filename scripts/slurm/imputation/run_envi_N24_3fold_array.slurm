#!/bin/bash
#SBATCH --job-name=envi_N24
#SBATCH --output=logs/envi_N24/%x_%A_%a.out
#SBATCH --error=logs/envi_N24/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-3%3

set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage:"
    echo "sbatch $0 DONOR EXPERIMENT EXPERIMENT_LABEL REFERENCE"
    exit 1
fi

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

SCRIPT="${PROJECT_ROOT}/src/imputation/ENVI/run_envi_fold.py"

DONOR="$1"
EXPERIMENT="$2"
EXPERIMENT_LABEL="$3"
REFERENCE="$4"

FOLD="${SLURM_ARRAY_TASK_ID}"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_full/gene_folds_200_100/${DONOR}"

OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_full/ENVI_N24/${DONOR}"

OBSERVED="${FOLD_DIR}/fold_${FOLD}_observed_genes.h5ad"
HELDOUT="${FOLD_DIR}/fold_${FOLD}_heldout_genes.h5ad"

mkdir -p \
    "${PROJECT_ROOT}/logs/envi_N24" \
    "${OUTPUT_ROOT}"

cd "${PROJECT_ROOT}"


# ============================================================
# INPUT VALIDATION
# ============================================================

if [[ ! -f "${OBSERVED}" ]]; then
    echo "ERROR: observed target missing:"
    echo "${OBSERVED}"
    exit 1
fi

if [[ ! -f "${HELDOUT}" ]]; then
    echo "ERROR: held-out truth missing:"
    echo "${HELDOUT}"
    exit 1
fi

if [[ ! -f "${REFERENCE}" ]]; then
    echo "ERROR: reference missing:"
    echo "${REFERENCE}"
    exit 1
fi


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
# ENVI ENVIRONMENT
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
    echo "ERROR: no working ENVI environment found."
    exit 1
fi


# ============================================================
# LOG CONFIGURATION
# ============================================================

echo "===================================================================================================="
echo "ENVI N24 THREE-FOLD OOF IMPUTATION"
echo "===================================================================================================="

echo "Job ID          : ${SLURM_JOB_ID}"
echo "Array job       : ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task      : ${SLURM_ARRAY_TASK_ID}"

echo
echo "Target donor    : ${DONOR}"
echo "Fold            : ${FOLD}"

echo
echo "Experiment      : ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"

echo
echo "Observed target : ${OBSERVED}"
echo "Held-out truth  : ${HELDOUT}"
echo "Reference       : ${REFERENCE}"
echo "Output root     : ${OUTPUT_ROOT}"

echo
echo "Python          : ${PYTHON}"
echo "Node            : $(hostname)"
echo "Started         : $(date)"

echo
echo "Training steps  : 16000"
echo "Batch size      : 128"
echo "num_HVG         : 2048"
echo "num_cov_genes   : 64"
echo "k_nearest       : 8"
echo "covet_batch_size: 256"
echo "spatial_dist    : pois"
echo "sc_dist         : nb"
echo "stable_eps      : 1e-6"

echo
echo "CPUs            : ${SLURM_CPUS_PER_TASK}"
echo "Memory          : 200G"
echo "Fold array      : 1-3%3"

echo "===================================================================================================="


# ============================================================
# RUN ENVI
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
echo "Donor     : ${DONOR}"
echo "Experiment: ${EXPERIMENT}"
echo "Fold      : ${FOLD}"
echo "Completed : $(date)"
echo "===================================================================================================="
