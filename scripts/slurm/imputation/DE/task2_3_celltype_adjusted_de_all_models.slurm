#!/bin/bash
#SBATCH --job-name=t23_ct_de
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/task2_3_celltype_de/%j.out
#SBATCH --error=logs/task2_3_celltype_de/%j.err

set -euo pipefail

cd /beegfs/labs/hulab/projects/mjabin/GeneBridge

RSCRIPT=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript

echo "================================================================================"
echo "TASK 2.3 — ALL 6 MODELS / N24 OOF / CELL-TYPE-ADJUSTED DE"
echo "================================================================================"
echo "Started: $(date)"
echo "Rscript: ${RSCRIPT}"
echo

"${RSCRIPT}" \
src/imputation/DE/task2_3_ntc_dysregulation/05_celltype_adjusted_de_all_models.R

echo
echo "================================================================================"
echo "SLURM STATUS: PASS"
echo "Finished: $(date)"
echo "================================================================================"
