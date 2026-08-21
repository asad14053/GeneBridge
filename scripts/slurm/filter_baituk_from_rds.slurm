#!/bin/bash
#SBATCH --job-name=filter_baituk_rds
#SBATCH --output=logs/filter_baituk_rds_%j.out
#SBATCH --error=logs/filter_baituk_rds_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"
SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/05_filter_baituk_from_rds.R"

export R_LIBS_USER="${PROJECT_ROOT}/software/R_libs"
export HDF5_USE_FILE_LOCKING=FALSE

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Create filtered Baituk H5AD from original RDS"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Rscript: ${RSCRIPT}"
echo "Script: ${SCRIPT}"
echo "============================================================"

"${RSCRIPT}" --version

"${RSCRIPT}" -e "
script <- '${SCRIPT}'
parse(file = script)
cat('R syntax check: PASSED\n')
"

"${RSCRIPT}" "${SCRIPT}"

echo ""
echo "============================================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "============================================================"
