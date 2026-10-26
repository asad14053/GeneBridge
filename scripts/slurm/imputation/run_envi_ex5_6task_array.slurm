#!/bin/bash
#SBATCH --job-name=envi_ex5
#SBATCH --output=logs/envi_N24/%x_%A_%a.out
#SBATCH --error=logs/envi_N24/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-6%6

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
SCRIPT="${PROJECT_ROOT}/src/imputation/ENVI/run_envi_fold.py"

FOLD_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/gene_folds_200_100"
REF_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/snrnaseq_target_references"
OUTPUT_ROOT_BASE="${PROJECT_ROOT}/outputs/imputation_full/ENVI_N24"

TASK="${SLURM_ARRAY_TASK_ID}"

# ============================================================
# Map 6 tasks to donor + fold
# ============================================================

case "${TASK}" in
    1)
        DONOR="Br6432"
        FOLD=1
        REFERENCE="${REF_ROOT}/Br6432/seq_data_huuki_snrna_Br6432_pseudocount.h5ad"
        ;;
    2)
        DONOR="Br6432"
        FOLD=2
        REFERENCE="${REF_ROOT}/Br6432/seq_data_huuki_snrna_Br6432_pseudocount.h5ad"
        ;;
    3)
        DONOR="Br6432"
        FOLD=3
        REFERENCE="${REF_ROOT}/Br6432/seq_data_huuki_snrna_Br6432_pseudocount.h5ad"
        ;;
    4)
        DONOR="Br8667"
        FOLD=1
        REFERENCE="${REF_ROOT}/Br8667/seq_data_huuki_snrna_Br8667.h5ad"
        ;;
    5)
        DONOR="Br8667"
        FOLD=2
        REFERENCE="${REF_ROOT}/Br8667/seq_data_huuki_snrna_Br8667.h5ad"
        ;;
    6)
        DONOR="Br8667"
        FOLD=3
        REFERENCE="${REF_ROOT}/Br8667/seq_data_huuki_snrna_Br8667.h5ad"
        ;;
    *)
        echo "ERROR: Invalid SLURM_ARRAY_TASK_ID=${TASK}"
        exit 1
        ;;
esac

EXPERIMENT="ex5"
EXPERIMENT_LABEL="Experiment 5"

FOLD_DIR="${FOLD_ROOT}/${DONOR}"
OUTPUT_ROOT="${OUTPUT_ROOT_BASE}/${DONOR}"

OBSERVED="${FOLD_DIR}/fold_${FOLD}_observed_genes.h5ad"
HELDOUT="${FOLD_DIR}/fold_${FOLD}_heldout_genes.h5ad"

mkdir -p "${OUTPUT_ROOT}" logs/envi_N24

# ============================================================
# Validate inputs
# ============================================================

for FILE in \
    "${REFERENCE}" \
    "${OBSERVED}" \
    "${HELDOUT}"
do
    if [[ ! -f "${FILE}" ]]; then
        echo "ERROR: Missing file:"
        echo "${FILE}"
        exit 1
    fi
done

cd "${PROJECT_ROOT}"

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
# Find working ENVI environment
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
    exit 1
fi

# ============================================================
# Log configuration
# ============================================================

echo "===================================================================================================="
echo "N24 ENVI 6-TASK EX-5 BENCHMARK"
echo "===================================================================================================="
echo "Job ID          : ${SLURM_JOB_ID}"
echo "Array task      : ${SLURM_ARRAY_TASK_ID}"
echo "Node            : $(hostname)"
echo "Started         : $(date)"
echo
echo "Target donor    : ${DONOR}"
echo "Experiment      : ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"
echo "Fold            : ${FOLD}"
echo
echo "Reference:"
echo "${REFERENCE}"
echo
echo "Observed target:"
echo "${OBSERVED}"
echo
echo "Held-out truth:"
echo "${HELDOUT}"
echo
echo "Output root:"
echo "${OUTPUT_ROOT}"
echo
echo "Python          : ${PYTHON}"
echo "Training steps  : 16000"
echo "Batch size      : 128"
echo "num_HVG         : 2048"
echo "num_cov_genes   : 64"
echo "k_nearest       : 8"
echo "covet_batch_size: 256"
echo "spatial_dist    : pois"
echo "sc_dist         : nb"
echo "stable_eps      : 1e-6"
echo "Seed            : $((8667 + TASK))"
echo
echo "CPUs            : ${SLURM_CPUS_PER_TASK}"
echo "Memory          : 200G"
echo "===================================================================================================="

# ============================================================
# Run one donor × one fold
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
    --seed "$((8667 + TASK))" \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --k-spatial 15 \
    --ssim-grid-size 128 \
    --n-clusters 10 \
    --n-pcs 30 \
    --plot-genes 10

echo
echo "----------------------------------------------------------------------------------------------------"
echo "FOLD STATUS: PASS"
echo "Donor     : ${DONOR}"
echo "Experiment: ${EXPERIMENT}"
echo "Fold      : ${FOLD}"
echo "Completed : $(date)"
echo "----------------------------------------------------------------------------------------------------"

echo
echo "FINAL STATUS: PASS"
