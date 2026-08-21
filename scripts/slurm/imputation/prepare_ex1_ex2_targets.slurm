#!/bin/bash
#SBATCH --job-name=prep_ex1_ex2
#SBATCH --output=logs/prep_ex1_ex2_%j.out
#SBATCH --error=logs/prep_ex1_ex2_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=160G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"
mkdir -p logs

echo "============================================================"
echo "Prepare Ex1 NTC + Ex2 SCZ production targets"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

"${PYTHON}" -u \
    src/imputation/production/00_prepare_ex1_ex2_targets.py

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
