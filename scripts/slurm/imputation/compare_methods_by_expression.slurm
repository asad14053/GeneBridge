#!/bin/bash
#SBATCH --job-name=compare_methods_expr
#SBATCH --output=logs/compare_methods_expr_%j.out
#SBATCH --error=logs/compare_methods_expr_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

PROJECT_ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

cd "$PROJECT_ROOT"

mkdir -p logs

echo "============================================================"
echo "Direct cross-method comparison by expression"
echo "Job: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================================"

"$PYTHON" -u \
  src/imputation/visualization/05_compare_methods_by_expression.py

echo
echo "============================================================"
echo "Completed: $(date)"
echo "============================================================"
