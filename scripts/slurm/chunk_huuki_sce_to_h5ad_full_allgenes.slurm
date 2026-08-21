#!/bin/bash
#SBATCH --job-name=huuki_sce_h5ad
#SBATCH --output=logs/huuki_sce_h5ad_%j.out
#SBATCH --error=logs/huuki_sce_h5ad_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

# ============================================================
# Paths
# ============================================================

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

ENV_ROOT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy"

RSCRIPT="${ENV_ROOT}/bin/Rscript"

SCRIPT="src/snrnaseq_branch/chunk_huuki_sce_to_h5ad_full_allgenes.R"


# ============================================================
# Job information
# ============================================================

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Working directory: ${PROJECT_ROOT}"
echo "Environment: ${ENV_ROOT}"
echo "========================================"

cd "${PROJECT_ROOT}"

mkdir -p logs


# ============================================================
# Check Rscript
# ============================================================

if [ ! -x "${RSCRIPT}" ]; then
    echo "ERROR: Rscript not found in genebridge_scanpy environment:"
    echo "${RSCRIPT}"
    exit 1
fi

echo ""
echo "Rscript:"
echo "${RSCRIPT}"

"${RSCRIPT}" --version


# ============================================================
# Run script
# ============================================================

echo ""
echo "========================================"
echo "Starting Huuki SCE -> H5AD full all-genes conversion"
echo "========================================"

"${RSCRIPT}" "${SCRIPT}" full 200 0


# ============================================================
# Done
# ============================================================

echo ""
echo "========================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "========================================"