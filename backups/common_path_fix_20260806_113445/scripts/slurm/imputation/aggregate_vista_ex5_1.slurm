#!/bin/bash
#SBATCH --job-name=vista_ex5_1_aggregate
#SBATCH --output=logs/vista_ex5_1_aggregate_%j.out
#SBATCH --error=logs/vista_ex5_1_aggregate_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"
OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

SCRIPT="${PROJECT_ROOT}/src/imputation/common/aggregate_three_folds.py"

# Support the earlier project location without modifying either file.
if [[ ! -f "${SCRIPT}" ]]; then
    SCRIPT="${PROJECT_ROOT}/src/imputation_beta/common/aggregate_three_folds.py"
fi

cd "${PROJECT_ROOT}"

mkdir -p \
  logs \
  "${OUTPUT_ROOT}/ex5_1/vista/combined"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

export PYTHONPATH="${PROJECT_ROOT}/src/imputation/common:${PROJECT_ROOT}/src/imputation_beta/common:${PYTHONPATH:-}"

echo "===================================================================================================="
echo "Aggregate VISTA Experiment 5.1 three-fold results"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 160G"
echo "Python: ${PYTHON}"
echo "Aggregation script: ${SCRIPT}"
echo "Fold directory: ${FOLD_DIR}"
echo "Output root: ${OUTPUT_ROOT}"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${SCRIPT}" || {
    echo "ERROR: aggregate_three_folds.py was not found in either expected location."
    exit 1
}

for FOLD in 1 2 3
do
    FOLD_OUTPUT="${OUTPUT_ROOT}/ex5_1/vista/fold_${FOLD}"

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

"${PYTHON}" --version
"${PYTHON}" -m py_compile "${SCRIPT}"

echo "Aggregation script syntax check: PASSED"

srun --cpu-bind=cores \
  "${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment ex5_1 \
    --experiment-label "Experiment 5.1" \
    --model vista \
    --model-label "VISTA" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --n-folds 3 \
    --n-clusters 10 \
    --n-pcs 30 \
    --seed 8667

echo
echo "===================================================================================================="
echo "VISTA Experiment 5.1 aggregation finished successfully"
echo "Completed: $(date)"
echo "===================================================================================================="
