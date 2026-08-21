#!/bin/bash
#SBATCH --job-name=prepare_ex5_1
#SBATCH --output=logs/prepare_ex5_1_%j.out
#SBATCH --error=logs/prepare_ex5_1_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/common/prepare_ex5_1_manifest.py"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export HDF5_USE_FILE_LOCKING=FALSE

cd "${PROJECT_ROOT}"

mkdir -p \
  logs \
  outputs/imputation_beta/Br8667/experiment_5_1/input_validation

echo "===================================================================================================="
echo "Experiment 5.1 preflight and manifest"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Memory: 120G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found: ${PYTHON}"
    exit 1
}

test -f "${SCRIPT}" || {
    echo "ERROR: Python script not found: ${SCRIPT}"
    exit 1
}

"${PYTHON}" --version
"${PYTHON}" -m py_compile "${SCRIPT}"

echo "Python syntax check: PASSED"

"${PYTHON}" -u "${SCRIPT}"

echo
echo "Completed: $(date)"
