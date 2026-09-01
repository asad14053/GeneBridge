#!/bin/bash
#SBATCH --job-name=ct_pba_de
#SBATCH --output=logs/celltype_de/11_pba_%j.out
#SBATCH --error=logs/celltype_de/11_pba_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00

set -euo pipefail

ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
RSCRIPT=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript

cd "$ROOT"

echo "============================================================"
echo "GLOBAL CELL-TYPE PAPER vs BEFORE vs AFTER DE"
echo "Job: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

"$RSCRIPT" \
src/imputation/DE/paper_concordance/11_celltype_global_paper_before_after_DE.R

echo
echo "============================================================"
echo "Finished successfully: $(date)"
echo "============================================================"
