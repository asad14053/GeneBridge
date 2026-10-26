#!/bin/bash
#SBATCH --job-name=aivc_audit
#SBATCH --output=logs/aivc_presplit_audit_%j.out
#SBATCH --error=logs/aivc_presplit_audit_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

mkdir -p logs
mkdir -p outputs/AIVC_beta/presplit_audit

echo "======================================================================"
echo "AIVC BETA PRE-SPLIT AUDIT"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "======================================================================"

"${PYTHON}" \
    src/AIVC/data_prep/00_presplit_audit.py

echo "======================================================================"
echo "Finished: $(date)"
echo "======================================================================"
