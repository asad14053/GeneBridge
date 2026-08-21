#!/bin/bash
#SBATCH --job-name=ex5_3_huuki9b
#SBATCH --output=logs/ex5_3_huuki9b_%j.out
#SBATCH --error=logs/ex5_3_huuki9b_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/05_3_create_huuki_9brains_excluding_xenium.py"

export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "===================================================================================================="
echo "Experiment 5.3: remove all Huuki variants associated with the Xenium brain"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "===================================================================================================="

"${PYTHON}" --version
"${PYTHON}" -m py_compile "${SCRIPT}"

echo "Python syntax check: PASSED"
echo
echo "Starting Experiment 5.3..."
echo

"${PYTHON}" -u "${SCRIPT}"

echo
echo "===================================================================================================="
echo "Experiment 5.3 completed successfully"
echo "Finished: $(date)"
echo "===================================================================================================="
