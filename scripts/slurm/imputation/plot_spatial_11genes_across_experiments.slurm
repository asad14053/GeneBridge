#!/bin/bash
#SBATCH --job-name=spatial11_exps
#SBATCH --output=logs/spatial11_exps_%A_%a.out
#SBATCH --error=logs/spatial11_exps_%A_%a.err
#SBATCH --array=1-6%6
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

PROJECT_ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

cd "$PROJECT_ROOT"

mkdir -p logs

METHODS=(
  vista
  gimvi
  tangram
  envi
  spage
  transimpspa
)

METHOD="${METHODS[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "============================================================"
echo "11-gene spatial panel across experiments"
echo "Method: ${METHOD}"
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

"$PYTHON" -u \
  src/imputation/visualization/02_spatial_11genes_across_experiments.py \
  --method "$METHOD"

echo
echo "============================================================"
echo "Completed: $(date)"
echo "Method: ${METHOD}"
echo "============================================================"
