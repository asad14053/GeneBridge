#!/bin/bash
#SBATCH --job-name=cell_layer_annotation
#SBATCH --output=logs/cell_layer_annotation_%j.out
#SBATCH --error=logs/cell_layer_annotation_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="/beegfs/labs/hulab/projects/mjabin/GeneBridge/src/xenium_branch/03_annotation/02_banksy_celltype_annotation.py"

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
