#!/bin/bash
#SBATCH --job-name=metric_expr
#SBATCH --output=logs/metric_expr_%A_%a.out
#SBATCH --error=logs/metric_expr_%A_%a.err
#SBATCH --array=1-6%6
#SBATCH --time=02:00:00
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
echo "Expression-stratified imputation metrics"
echo "Method: ${METHOD}"
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Started: $(date)"
echo "============================================================"

"$PYTHON" -u \
  src/imputation/visualization/04_metric_distributions_by_expression.py \
  --method "$METHOD"

echo
echo "============================================================"
echo "Completed: $(date)"
echo "Method: ${METHOD}"
echo "============================================================"
