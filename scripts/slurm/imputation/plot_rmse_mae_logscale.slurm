#!/bin/bash
#SBATCH --job-name=rmse_mae_log
#SBATCH --output=logs/rmse_mae_log_%j.out
#SBATCH --error=logs/rmse_mae_log_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

PROJECT_ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

cd "$PROJECT_ROOT"

mkdir -p logs

echo "============================================================"
echo "RMSE / MAE log-scale cross-method comparison"
echo "Job: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================================"

"$PYTHON" -u \
  src/imputation/visualization/06_rmse_mae_logscale_comparison.py

echo
echo "============================================================"
echo "Completed: $(date)"
echo "============================================================"
