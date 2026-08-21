#!/bin/bash
#SBATCH --job-name=baituk_huuki_filter
#SBATCH --output=logs/baituk_huuki_filter_%j.out
#SBATCH --error=logs/baituk_huuki_filter_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/04_create_overlap_filtered_h5ads.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Create filtered Baituk and Huuki overlap H5AD files"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "============================================================"

if [ ! -x "${PYTHON}" ]; then
    echo "ERROR: Python not found: ${PYTHON}"
    exit 1
fi

if [ ! -f "${SCRIPT}" ]; then
    echo "ERROR: Python script not found: ${SCRIPT}"
    exit 1
fi

"${PYTHON}" --version

"${PYTHON}" -c "
import sys
import anndata
import pandas
import numpy

print('Python executable:', sys.executable)
print('anndata:', anndata.__version__)
print('pandas:', pandas.__version__)
print('numpy:', numpy.__version__)
"

echo ""
echo "Starting overlap filtering..."

"${PYTHON}" -u "${SCRIPT}"

echo ""
echo "============================================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "============================================================"
