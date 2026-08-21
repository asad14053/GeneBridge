#!/bin/bash
#SBATCH --job-name=batiuk_huuki_genes
#SBATCH --output=logs/batiuk_huuki_genes_%j.out
#SBATCH --error=logs/batiuk_huuki_genes_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/02_compare_batiuk_huuki_genes.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Compare Batiuk and Huuki genes"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Project root: ${PROJECT_ROOT}"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "============================================================"

if [ ! -x "${PYTHON}" ]; then
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
fi

if [ ! -f "${SCRIPT}" ]; then
    echo "ERROR: Python script not found:"
    echo "${SCRIPT}"
    exit 1
fi

"${PYTHON}" --version
"${PYTHON}" -c "import sys, anndata, pandas; print('Python:', sys.executable); print('anndata:', anndata.__version__); print('pandas:', pandas.__version__)"

echo ""
echo "Starting Batiuk–Huuki gene comparison..."
echo ""

"${PYTHON}" -u "${SCRIPT}"

echo ""
echo "============================================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "============================================================"
