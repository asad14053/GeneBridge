#!/bin/bash
#SBATCH --job-name=inspect_xenium_Br8667
#SBATCH --output=logs/inspect_xenium_Br8667_%j.out
#SBATCH --error=logs/inspect_xenium_Br8667_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="src/imputation/05_inspect_xenium_Br8667.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "========================================"

"${PYTHON}" --version

"${PYTHON}" -u "${SCRIPT}"

echo "========================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "========================================"
