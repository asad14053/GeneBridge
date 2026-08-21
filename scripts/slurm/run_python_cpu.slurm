#!/bin/bash
#SBATCH --job-name=merge_huuki
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch $0 <python_script> [script_arguments...]"
    exit 1
fi

SCRIPT="$1"
shift

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

echo ""
echo "Starting Python script..."
echo ""

"${PYTHON}" -u "${SCRIPT}" "$@"

echo ""
echo "========================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "========================================"