#!/bin/bash
#SBATCH --job-name=check_ref_align
#SBATCH --output=logs/check_ref_align_%j.out
#SBATCH --error=logs/check_ref_align_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/07_check_three_reference_target_alignment.py"

export HDF5_USE_FILE_LOCKING=FALSE

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Check three references against one Br8667 Xenium target"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "============================================================"

"${PYTHON}" -m py_compile "${SCRIPT}"
echo "Python syntax check: PASSED"

"${PYTHON}" -u "${SCRIPT}"

echo
echo "Completed successfully: $(date)"
