#!/bin/bash
#SBATCH --job-name=vista_umap
#SBATCH --output=logs/vista_umap_%j.out
#SBATCH --error=logs/vista_umap_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

PROJECT_ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

cd "$PROJECT_ROOT"

mkdir -p logs

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "============================================================"
echo "VISTA shared-basis UMAP"
echo "Job: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================================"

"$PYTHON" -u \
  src/imputation/visualization/09_umap_vista_across_experiments.py

echo
echo "============================================================"
echo "Completed: $(date)"
echo "============================================================"
