#!/bin/bash
#SBATCH --job-name=imputation_agg_v2
#SBATCH --output=logs/imputation_agg_v2_%j.out
#SBATCH --error=logs/imputation_agg_v2_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

AGGREGATOR="${PROJECT_ROOT}/src/imputation/common/aggregate_three_folds_v2.py"

EXPERIMENT="${1:?Missing experiment key, for example ex5}"
EXPERIMENT_LABEL="${2:?Missing experiment label}"
MODEL="${3:?Missing model key, for example vista}"
MODEL_LABEL="${4:?Missing model label}"
COMBINED_DIR_NAME="${5:-combined_v2}"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"
OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

export PYTHONPATH="${PROJECT_ROOT}/src/imputation/common:${PYTHONPATH:-}"

cd "${PROJECT_ROOT}"

mkdir -p logs

echo "===================================================================================================="
echo "Shared three-fold imputation aggregation v2"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo
echo "Memory: 160G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Experiment: ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"
echo "Model: ${MODEL}"
echo "Model label: ${MODEL_LABEL}"
echo "Combined directory: ${COMBINED_DIR_NAME}"
echo "Aggregator: ${AGGREGATOR}"
echo "Common directory: ${PROJECT_ROOT}/src/imputation/common"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${AGGREGATOR}" || {
    echo "ERROR: Shared aggregator not found:"
    echo "${AGGREGATOR}"
    exit 1
}

for FOLD in 1 2 3
do
    FOLD_OUTPUT="${OUTPUT_ROOT}/${EXPERIMENT}/${MODEL}/fold_${FOLD}"

    for REQUIRED_FILE in \
        "${FOLD_OUTPUT}/run_complete.flag" \
        "${FOLD_OUTPUT}/predicted_heldout_genes.h5ad" \
        "${FOLD_OUTPUT}/gene_level_metrics.csv" \
        "${FOLD_OUTPUT}/fold_level_metrics.csv"
    do
        test -f "${REQUIRED_FILE}" || {
            echo "ERROR: Required fold output not found:"
            echo "${REQUIRED_FILE}"
            exit 1
        }
    done
done

"${PYTHON}" -m py_compile \
  "${PROJECT_ROOT}/src/imputation/common/benchmark_evaluation_v2.py" \
  "${AGGREGATOR}"

echo "Shared module syntax check: PASSED"

srun --cpu-bind=cores \
  "${PYTHON}" -u "${AGGREGATOR}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --experiment-label "${EXPERIMENT_LABEL}" \
    --model "${MODEL}" \
    --model-label "${MODEL_LABEL}" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --n-folds 3 \
    --n-clusters 10 \
    --n-pcs 30 \
    --seed 8667 \
    --combined-dir-name "${COMBINED_DIR_NAME}"

echo
echo "Completed successfully: $(date)"
