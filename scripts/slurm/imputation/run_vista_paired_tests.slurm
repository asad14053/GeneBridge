#!/bin/bash
#SBATCH --job-name=vista_paired_tests
#SBATCH --output=logs/vista_paired_tests_%j.out
#SBATCH --error=logs/vista_paired_tests_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/common/paired_tests_vista_experiments.py"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

cd "${PROJECT_ROOT}"

mkdir -p \
  logs \
  outputs/imputation_beta/Br8667/statistical_tests/vista_experiments

echo "===================================================================================================="
echo "VISTA paired t-tests and 5% TOST equivalence tests"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Memory: 160G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${SCRIPT}" || {
    echo "ERROR: Statistical test script not found:"
    echo "${SCRIPT}"
    exit 1
}

for EXPERIMENT_DIR in \
    "${PROJECT_ROOT}/outputs/imputation_beta/Br8667/ex5/vista" \
    "${PROJECT_ROOT}/outputs/imputation_beta/Br8667/ex5_1/vista" \
    "${PROJECT_ROOT}/outputs/imputation_beta/Br8667/ex5_3/vista"
do
    test -d "${EXPERIMENT_DIR}" || {
        echo "ERROR: Experiment output directory not found:"
        echo "${EXPERIMENT_DIR}"
        exit 1
    }
done

"${PYTHON}" --version

"${PYTHON}" -m py_compile "${SCRIPT}"

echo "Python syntax check: PASSED"

srun --cpu-bind=cores \
  "${PYTHON}" -u "${SCRIPT}"

echo
echo "Completed: $(date)"
