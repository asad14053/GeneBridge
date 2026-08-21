#!/bin/bash
#SBATCH --job-name=tangram_3fold
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-3%3

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage:"
    echo "sbatch --job-name=tangram_ex5 $0 EXPERIMENT EXPERIMENT_LABEL REFERENCE"
    exit 1
fi

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/Tangram/run_tangram_fold.py"
FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"
OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

EXPERIMENT="$1"
EXPERIMENT_LABEL="$2"
REFERENCE="$3"

FOLD="${SLURM_ARRAY_TASK_ID}"

cd "${PROJECT_ROOT}"

mkdir -p logs "${OUTPUT_ROOT}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

echo "===================================================================================================="
echo "Tangram three-fold benchmark"
echo "===================================================================================================="
echo "Job: ${SLURM_JOB_ID}"
echo "Array job: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Task/fold: ${FOLD}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo
echo "Experiment: ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"
echo "Reference: ${REFERENCE}"
echo "Tangram mode: clusters"
echo "Epochs: 500"
echo "Device: cpu"
echo "Memory: 200G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Array policy: 1-3%3"
echo "===================================================================================================="

"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --experiment-label "${EXPERIMENT_LABEL}" \
    --fold "${FOLD}" \
    --reference "${REFERENCE}" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --tangram-mode clusters \
    --cluster-label auto \
    --min-reference-cells-per-cluster 10 \
    --density-prior uniform \
    --normalization-target-sum 10000 \
    --run-mode full \
    --full-epochs 500 \
    --device cpu \
    --seed 8667 \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --k-spatial 15 \
    --ssim-grid-size 128 \
    --n-clusters 10 \
    --n-pcs 30 \
    --plot-genes 10 \
    --skip-map-save

echo
echo "===================================================================================================="
echo "SUCCESS"
echo "Experiment: ${EXPERIMENT}"
echo "Fold: ${FOLD}"
echo "Completed: $(date)"
echo "===================================================================================================="
