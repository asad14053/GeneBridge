#!/bin/bash
#SBATCH --job-name=prep_ex1_ex2_fast
#SBATCH --output=logs/prep_ex1_ex2_fast_%j.out
#SBATCH --error=logs/prep_ex1_ex2_fast_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="src/imputation/production/00_prepare_ex1_ex2_targets_fast.py"

cd "${ROOT}"
mkdir -p logs

echo "============================================================"
echo "Ex1 NTC + Ex2 SCZ target preparation"
echo "RAM: 200G"
echo "CPUs: 16"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

"${PYTHON}" -u "${SCRIPT}"

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
