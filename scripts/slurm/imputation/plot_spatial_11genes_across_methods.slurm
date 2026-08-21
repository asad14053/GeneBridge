#!/bin/bash
#SBATCH --job-name=spatial11_methods
#SBATCH --output=logs/spatial11_methods_%A_%a.out
#SBATCH --error=logs/spatial11_methods_%A_%a.err
#SBATCH --array=1-3%3
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

PROJECT_ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

cd "$PROJECT_ROOT"

mkdir -p logs

EXPERIMENTS=(
  ex5
  ex5_1
  ex5_3
)

EXPERIMENT="${EXPERIMENTS[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "============================================================"
echo "11-gene spatial panel across methods"
echo "Experiment: ${EXPERIMENT}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

"$PYTHON" -u \
  src/imputation/visualization/01_spatial_11genes_across_methods.py \
  --experiment "$EXPERIMENT"

echo
echo "============================================================"
echo "Completed: $(date)"
echo "Experiment: ${EXPERIMENT}"
echo "============================================================"
