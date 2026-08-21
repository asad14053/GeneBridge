#!/bin/bash
#SBATCH --job-name=huuki_snrna_meta
#SBATCH --output=logs/huuki_snrna_meta_%j.out
#SBATCH --error=logs/huuki_snrna_meta_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

# ============================================================
# Paths
# ============================================================

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"

SCRIPT="src/snrnaseq_branch/02_metadata_qc/01_export_huuki_snrna_metadata.R"


# ============================================================
# Job information
# ============================================================

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Working directory: ${PROJECT_ROOT}"
echo "========================================"

cd "${PROJECT_ROOT}"

mkdir -p logs


# ============================================================
# Environment check
# ============================================================

echo ""
echo "Rscript executable:"
echo "${RSCRIPT}"

echo ""
echo "R version:"
"${RSCRIPT}" --version


# ============================================================
# Run script
# ============================================================

echo ""
echo "========================================"
echo "Starting Huuki snRNA-seq metadata export"
echo "========================================"

"${RSCRIPT}" "${SCRIPT}"


# ============================================================
# Done
# ============================================================

echo ""
echo "========================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "========================================"