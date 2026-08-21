#!/bin/bash
#SBATCH --job-name=xenium_qc
#SBATCH --output=logs/xenium_qc_%j.out
#SBATCH --error=logs/xenium_qc_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

# ============================================================
# Paths
# ============================================================

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

SCRIPT="src/xenium_branch/02_qc/01_xenium_qc.py"


# ============================================================
# Job information
# ============================================================

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Working directory: ${PROJECT_ROOT}"
echo "========================================"


# ============================================================
# Move to GeneBridge
# ============================================================

cd "${PROJECT_ROOT}"

mkdir -p logs


# ============================================================
# Environment check
# ============================================================

echo ""
echo "Python executable:"
echo "${PYTHON}"

echo ""
echo "Python version:"
"${PYTHON}" --version

echo ""
echo "Python environment:"
"${PYTHON}" -c "import sys; print(sys.executable)"


# ============================================================
# Run Xenium QC
# ============================================================

echo ""
echo "========================================"
echo "Starting Xenium QC"
echo "========================================"

"${PYTHON}" -u "${SCRIPT}"


# ============================================================
# Done
# ============================================================

echo ""
echo "========================================"
echo "Xenium QC completed successfully"
echo "Completed: $(date)"
echo "========================================"