#!/bin/bash
#SBATCH --job-name=extract_baituk_genes
#SBATCH --output=logs/extract_baituk_genes_%j.out
#SBATCH --error=logs/extract_baituk_genes_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"
SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/01_extract_batiuk_genes.R"

export R_LIBS_USER="${PROJECT_ROOT}/software/R_libs"
export HDF5_USE_FILE_LOCKING=FALSE

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Extract Baituk genes"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Project root: ${PROJECT_ROOT}"
echo "Rscript: ${RSCRIPT}"
echo "Script: ${SCRIPT}"
echo "============================================================"

if [ ! -x "${RSCRIPT}" ]; then
    echo "ERROR: Rscript executable not found:"
    echo "${RSCRIPT}"
    exit 1
fi

if [ ! -f "${SCRIPT}" ]; then
    echo "ERROR: R script not found:"
    echo "${SCRIPT}"
    exit 1
fi

"${RSCRIPT}" --version

echo ""
echo "Checking R script syntax..."

"${RSCRIPT}" -e "
parse(file='${SCRIPT}')
cat('R syntax check: PASSED\n')
"

echo ""
echo "Starting Baituk gene extraction..."
echo ""

"${RSCRIPT}" "${SCRIPT}"

echo ""
echo "============================================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "============================================================"
