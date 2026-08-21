#!/bin/bash
#SBATCH --job-name=transimpspa_3fold
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-3%3

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage:"
    echo "sbatch $0 EXPERIMENT EXPERIMENT_LABEL REFERENCE"
    exit 1
fi

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/TransImp/run_transimpspa_fold.py"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"

OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

EXPERIMENT="$1"
EXPERIMENT_LABEL="$2"
REFERENCE="$3"

FOLD="${SLURM_ARRAY_TASK_ID}"

cd "${PROJECT_ROOT}"

mkdir -p logs

export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export PYTHONPATH="${PROJECT_ROOT}/src/imputation/TransImp/tranSpa:${PROJECT_ROOT}/src/imputation/common:${PYTHONPATH:-}"

echo "===================================================================================================="
echo "TransImpSpa three-fold benchmark"
echo "===================================================================================================="
echo "Job: ${SLURM_JOB_ID}"
echo "Fold: ${FOLD}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo
echo "Experiment: ${EXPERIMENT}"
echo "Reference: ${REFERENCE}"
echo
echo "signature_mode: cell"
echo "mapping_mode: lowrank"
echo "mapping_lowdim: 256"
echo "wt_spa: 1.0"
echo "spatial kNN: 6"
echo "epochs: 2000"
echo "lr: 0.01"
echo "weight_decay: 0.01"
echo "clip_max: 10"
echo
echo "Memory: 200G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "===================================================================================================="

"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --experiment-label "${EXPERIMENT_LABEL}" \
    --fold "${FOLD}" \
    --reference "${REFERENCE}" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --mapping-lowdim 256 \
    --epochs 2000 \
    --lr 0.01 \
    --weight-decay 0.01 \
    --clip-max 10 \
    --wt-spa 1.0 \
    --spatial-neighbors 6 \
    --seed 8667 \
    --cpus "${SLURM_CPUS_PER_TASK}"

echo
echo "===================================================================================================="
echo "SUCCESS"
echo "Model: TransImpSpa"
echo "Experiment: ${EXPERIMENT}"
echo "Fold: ${FOLD}"
echo "Completed: $(date)"
echo "===================================================================================================="
