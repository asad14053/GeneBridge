#!/bin/bash
#SBATCH --job-name=vista_ex5_3_aggregate
#SBATCH --output=logs/vista_ex5_3_aggregate_%j.out
#SBATCH --error=logs/vista_ex5_3_aggregate_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

AGGREGATOR="${PROJECT_ROOT}/src/imputation/common/aggregate_three_folds.py"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"

OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

cd "${PROJECT_ROOT}"

mkdir -p \
  logs \
  "${OUTPUT_ROOT}/ex5_3/vista/combined"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

export PYTHONPATH="${PROJECT_ROOT}/src/imputation/common:${PYTHONPATH:-}"

echo "===================================================================================================="
echo "Aggregate VISTA Experiment 5.3 three-fold results"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo
echo "Memory: 160G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Python: ${PYTHON}"
echo "Aggregator: ${AGGREGATOR}"
echo "Output root: ${OUTPUT_ROOT}"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${AGGREGATOR}" || {
    echo "ERROR: Aggregation script not found:"
    echo "${AGGREGATOR}"
    exit 1
}

for FOLD in 1 2 3
do
    FOLD_OUTPUT="${OUTPUT_ROOT}/ex5_3/vista/fold_${FOLD}"

    test -d "${FOLD_OUTPUT}" || {
        echo "ERROR: Fold output directory not found:"
        echo "${FOLD_OUTPUT}"
        exit 1
    }

    test -f "${FOLD_OUTPUT}/run_complete.flag" || {
        echo "ERROR: Fold ${FOLD} is not marked complete:"
        echo "${FOLD_OUTPUT}/run_complete.flag"
        exit 1
    }

    echo
    echo "Fold ${FOLD} output files:"
    find "${FOLD_OUTPUT}" \
      -maxdepth 1 \
      -type f \
      -printf "%f\n" \
      | sort
done

"${PYTHON}" --version
"${PYTHON}" -m py_compile "${AGGREGATOR}"

echo
echo "Aggregation script syntax check: PASSED"

srun --cpu-bind=cores \
  "${PYTHON}" -u "${AGGREGATOR}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment ex5_3 \
    --experiment-label "Experiment 5.3" \
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
echo "VISTA Experiment 5.3 aggregation completed successfully"
echo "Completed: $(date)"
echo "===================================================================================================="
