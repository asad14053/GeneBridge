#!/bin/bash
#SBATCH --job-name=vista_ex5_aggregate
#SBATCH --output=logs/vista_ex5_aggregate_%j.out
#SBATCH --error=logs/vista_ex5_aggregate_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/common/aggregate_three_folds.py"
FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"
OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

cd "${PROJECT_ROOT}"
mkdir -p logs

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MPLBACKEND=Agg

"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment ex5 \
    --experiment-label "Experiment 5" \
    --model vista \
    --model-label "VISTA" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --n-folds 3 \
    --n-clusters 10 \
    --n-pcs 30 \
    --seed 8667

echo "Finished: $(date)"
