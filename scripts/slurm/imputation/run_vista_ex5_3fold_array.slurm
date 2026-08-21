#!/bin/bash
#SBATCH --job-name=vista_ex5_3fold
#SBATCH --output=logs/vista_ex5_%A_%a.out
#SBATCH --error=logs/vista_ex5_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --array=1-3%3

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
export PYTHONPATH="${PROJECT_ROOT}/src/imputation/VISTA:${PROJECT_ROOT}/src/imputation/common:${PYTHONPATH:-}"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/VISTA/run_vista_fold.py"
REFERENCE="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/seq_data_huuki_snrna_Br8667_vista.h5ad"
FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"
OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

FOLD="${SLURM_ARRAY_TASK_ID}"

cd "${PROJECT_ROOT}"
mkdir -p logs "${OUTPUT_ROOT}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MPLBACKEND=Agg

echo "===================================================================================================="
echo "VISTA Experiment 5, fold ${FOLD}"
echo "Job: ${SLURM_JOB_ID}"
echo "Array job: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Reference: ${REFERENCE}"
echo "Fold directory: ${FOLD_DIR}"
echo "===================================================================================================="

srun --cpu-bind=cores \
    "${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment ex5 \
    --experiment-label "Experiment 5" \
    --fold "${FOLD}" \
    --reference "${REFERENCE}" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --reference-mode benchmark_300 \
    --run-mode full \
    --full-epochs 200 \
    --batch-size 128 \
    --n-latent 32 \
    --neighbor-size 20 \
    --seed 8667 \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --k-spatial 15 \
    --ssim-grid-size 128 \
    --n-clusters 10 \
    --n-pcs 30 \
    --plot-genes 10

echo "Finished: $(date)"
