#!/bin/bash
#SBATCH --job-name=ct_pb
#SBATCH --output=logs/celltype_pseudobulk/%A_%a.out
#SBATCH --error=logs/celltype_pseudobulk/%A_%a.err
#SBATCH --array=1-23%23
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00

set -euo pipefail

ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

cd "$ROOT"

echo "========================================"
echo "Cell-type pseudobulk"
echo "SLURM_JOB_ID:       ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID:${SLURM_ARRAY_TASK_ID}"
echo "Node:               $(hostname)"
echo "Started:            $(date)"
echo "========================================"

"$PYTHON" \
    src/imputation/DE/paper_concordance/10a_build_envi_celltype_pseudobulk_donor.py \
    --index "${SLURM_ARRAY_TASK_ID}"

echo
echo "========================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "========================================"
