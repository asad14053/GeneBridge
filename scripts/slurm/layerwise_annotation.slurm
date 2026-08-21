#!/bin/bash
#SBATCH --job-name=layer_annotation
#SBATCH --output=logs/layer_annotation_%j.out
#SBATCH --error=logs/layer_annotation_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="src/xenium_branch/03_annotation/01_layerwise_annotation.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "========================================"

"${PYTHON}" -u "${SCRIPT}" \
    --compression lzf \
    --max-cells-per-sample 0

echo "========================================"
echo "Completed: $(date)"
echo "========================================"
