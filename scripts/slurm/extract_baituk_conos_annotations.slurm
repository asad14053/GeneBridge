#!/bin/bash
#SBATCH --job-name=baituk_annotations
#SBATCH --output=logs/baituk_annotations_%j.out
#SBATCH --error=logs/baituk_annotations_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"
SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/04_baituk_celltype_annotation/01_extract_conos_annotations.R"

export R_LIBS_USER="${PROJECT_ROOT}/software/R_libs"
export HDF5_USE_FILE_LOCKING=FALSE

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Extract Baituk published cell annotations"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Rscript: ${RSCRIPT}"
echo "Script: ${SCRIPT}"
echo "============================================================"

"${RSCRIPT}" --version

"${RSCRIPT}" -e '
if (!requireNamespace("conos", quietly = TRUE)) {
    stop("The R package conos is not installed in genebridge_r.")
}
cat("conos version: ", as.character(packageVersion("conos")), "\n", sep="")
'

"${RSCRIPT}" -e "
parse(file='${SCRIPT}')
cat('R syntax check: PASSED\n')
"

"${RSCRIPT}" "${SCRIPT}"

echo
echo "Completed successfully: $(date)"
