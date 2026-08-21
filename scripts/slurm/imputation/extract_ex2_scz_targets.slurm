#!/bin/bash
#SBATCH --job-name=ex2_scz_targets
#SBATCH --output=logs/ex2_scz_targets_%j.out
#SBATCH --error=logs/ex2_scz_targets_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"
mkdir -p logs

echo "============================================================"
echo "Experiment 2 — Extract SCZ Xenium targets"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

"${PYTHON}" -u \
    src/imputation/production/01_extract_ex2_scz_targets.py

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
