#!/bin/bash
#SBATCH --job-name=huuki_visium_Br8667
#SBATCH --output=logs/huuki_visium_Br8667_%j.out
#SBATCH --error=logs/huuki_visium_Br8667_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/02_make_huuki_visium_Br8667_h5ad.py"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Huuki Visium Br8667 H5AD"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Project root: ${PROJECT_ROOT}"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "============================================================"

if [ ! -d "${PROJECT_ROOT}" ]; then
    echo "ERROR: Project directory not found:"
    echo "${PROJECT_ROOT}"
    exit 1
fi

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
"${PYTHON}" -c "import sys; print('Python executable:', sys.executable)"

echo ""
echo "Starting script..."

"${PYTHON}" -u "${SCRIPT}"

echo ""
echo "============================================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "============================================================"
