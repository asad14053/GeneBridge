#!/bin/bash
#SBATCH --job-name=ex5_3_huuki9
#SBATCH --output=logs/ex5_3_huuki9_%j.out
#SBATCH --error=logs/ex5_3_huuki9_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/05_3_create_huuki_9samples_excluding_xenium.py"

export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "===================================================================================================="
echo "Experiment 5.3: Huuki 10 samples minus Xenium-overlapping sample"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Project root: ${PROJECT_ROOT}"
echo "Python: ${PYTHON}"
echo "Script: ${SCRIPT}"
echo "===================================================================================================="

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

"${PYTHON}" -c "
import sys
import anndata
import pandas
import numpy
import scipy
import h5py

print('Python executable:', sys.executable)
print('anndata:', anndata.__version__)
print('pandas:', pandas.__version__)
print('numpy:', numpy.__version__)
print('scipy:', scipy.__version__)
print('h5py:', h5py.__version__)
"

"${PYTHON}" -m py_compile "${SCRIPT}"
echo "Python syntax check: PASSED"

echo
echo "Starting Experiment 5.3 data creation..."
echo

"${PYTHON}" -u "${SCRIPT}"

echo
echo "===================================================================================================="
echo "Experiment 5.3 completed successfully"
echo "Finished: $(date)"
echo "===================================================================================================="
