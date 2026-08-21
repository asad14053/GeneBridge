#!/bin/bash
#SBATCH --job-name=filter_huuki
#SBATCH --output=logs/filter_huuki_%j.out
#SBATCH --error=logs/filter_huuki_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/06_filter_huuki_h5ad.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Filter full Huuki H5AD"
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
