#!/bin/bash
#SBATCH --job-name=audit_ex5_matrix
#SBATCH --output=logs/audit_ex5_matrix_%j.out
#SBATCH --error=logs/audit_ex5_matrix_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/08_audit_expression_layers_experiments_5.py"

export HDF5_USE_FILE_LOCKING=FALSE

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"

"${PYTHON}" -m py_compile "${SCRIPT}"
"${PYTHON}" -u "${SCRIPT}"

echo "Completed: $(date)"
