#!/bin/bash
#SBATCH --job-name=filter_baituk
#SBATCH --output=logs/filter_baituk_%j.out
#SBATCH --error=logs/filter_baituk_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/05_filter_baituk_h5ad.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Filter existing full Baituk H5AD"
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

echo ""
echo "Completed successfully: $(date)"
