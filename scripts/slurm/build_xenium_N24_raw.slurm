#!/bin/bash
#SBATCH --job-name=xenium_N24_raw
#SBATCH --output=logs/xenium_N24_raw_%j.out
#SBATCH --error=logs/xenium_N24_raw_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

# ============================================================
# Paths
# ============================================================

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

SCRIPT="src/xenium_branch/01_build_datastructure/01_build_xenium_N24_raw.py"


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
# Move to project directory
# ============================================================

cd "${PROJECT_ROOT}"

mkdir -p logs


# ============================================================
# Check Python environment
# ============================================================

echo ""
echo "Python executable:"
echo "${PYTHON}"

echo ""
echo "Python version:"
"${PYTHON}" --version

echo ""
echo "Python environment check:"
"${PYTHON}" -c "import sys; print(sys.executable)"


# ============================================================
# Run Xenium N24 data structure build
# ============================================================

echo ""
echo "========================================"
echo "Starting Xenium N24 raw data structure build"
echo "========================================"

"${PYTHON}" -u "${SCRIPT}"


# ============================================================
# Done
# ============================================================

echo ""
echo "========================================"
echo "Xenium N24 build completed successfully"
echo "Completed: $(date)"
echo "========================================"