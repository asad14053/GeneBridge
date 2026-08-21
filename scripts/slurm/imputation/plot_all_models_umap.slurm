#!/bin/bash
#SBATCH --job-name=all_umap
#SBATCH --output=logs/all_umap_%A_%a.out
#SBATCH --error=logs/all_umap_%A_%a.err
#SBATCH --array=1-6%6
#SBATCH --time=06:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

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

METHOD="${METHODS[$((SLURM_ARRAY_TASK_ID - 1))]}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "============================================================"
echo "Shared-basis UMAP"
echo "Method: ${METHOD}"
echo "Array job: ${SLURM_ARRAY_JOB_ID}"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo "============================================================"

"$PYTHON" -u \
  src/imputation/visualization/10_umap_method_across_experiments.py \
  --method "$METHOD"

echo
echo "============================================================"
echo "Method complete: ${METHOD}"
echo "Completed: $(date)"
echo "============================================================"
