#!/bin/bash
#SBATCH --job-name=br8667_snrna
#SBATCH --output=logs/br8667_snrna_%j.out
#SBATCH --error=logs/br8667_snrna_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="src/snrnaseq_branch/create_brain_snrna_from_huuki_h5ad.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "========================================"

"${PYTHON}" --version

"${PYTHON}" -u "${SCRIPT}" \
    --brain-id Br8667 \
    --top-genes full \
    --compression lzf

echo "========================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "========================================"
