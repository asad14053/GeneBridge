#!/bin/bash
#SBATCH --job-name=t2.3_all300_QC
#SBATCH --output=logs/task2_3_de_all300_filterQC/%j.out
#SBATCH --error=logs/task2_3_de_all300_filterQC/%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00

set -euo pipefail

ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge
R_ENV=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r

cd "$ROOT"

echo "============================================================"
echo "TASK 2.3 — ALL 300 GENES + filterByExpr QC"
echo "============================================================"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo
echo "DE unit     : donor × spatial domain"
echo "Donors      : 24 (12 NTC + 12 SCZ)"
echo "Domains     : 7"
echo "Pseudobulks : 168"
echo "Genes       : all 300 OOF genes retained"
echo "Filter      : filterByExpr QC ONLY — no genes removed"
echo "Covariates  : predictions_smooth + Age + Sex + slide_id"
echo "Block       : BrNum"
echo "run_date    : NOT included"
echo "============================================================"
echo

"$R_ENV/bin/Rscript" \
    src/imputation/DE/task2_3_ntc_dysregulation/02_layer_adjusted_de_all_models_all300_filterQC.R

echo
echo "============================================================"
echo "SLURM STATUS: PASS"
echo "Finished: $(date)"
echo "============================================================"
