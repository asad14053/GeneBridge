#!/bin/bash
#SBATCH --job-name=visium_layers_N24
#SBATCH --output=logs/visium_layers_N24_%j.out
#SBATCH --error=logs/visium_layers_N24_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"
SCRIPT="src/visium_branch/03_annotation/01_visualize_matched_N24_visium_layers.R"

export R_LIBS_USER="${PROJECT_ROOT}/software/R_libs"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Rscript: ${RSCRIPT}"
echo "Script: ${SCRIPT}"
echo "========================================"

"${RSCRIPT}" --version

"${RSCRIPT}" "${SCRIPT}"

echo "========================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "========================================"
