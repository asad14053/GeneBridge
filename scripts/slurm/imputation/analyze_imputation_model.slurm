#!/bin/bash
#SBATCH --job-name=imputation_analysis
#SBATCH --output=logs/imputation_analysis_%j.out
#SBATCH --error=logs/imputation_analysis_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage:"
    echo "  sbatch $0 MODEL MODEL_LABEL [COMBINED_DIR]"
    echo
    echo "Example:"
    echo '  sbatch $0 gimvi "gimVI" combined_v2'
    exit 1
fi

MODEL="$1"
MODEL_LABEL="$2"
COMBINED_DIR="${3:-combined_v2}"

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

cd "${PROJECT_ROOT}"

mkdir -p logs

export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "===================================================================================================="
echo "STANDARDIZED IMPUTATION MODEL ANALYSIS"
echo "===================================================================================================="
echo "Model key:     ${MODEL}"
echo "Model label:   ${MODEL_LABEL}"
echo "Combined dir:  ${COMBINED_DIR}"
echo "Started:       $(date)"
echo "Host:          $(hostname)"
echo

echo "===================================================================================================="
echo "1. PAIRED T-TEST + TOST"
echo "===================================================================================================="

"${PYTHON}" -u \
  src/imputation/common/paired_tests_imputation_experiments.py \
  --model "${MODEL}" \
  --model-label "${MODEL_LABEL}" \
  --combined-dir "${COMBINED_DIR}"

echo
echo "===================================================================================================="
echo "2. FIXED ADVISOR 11-GENE COMPARISON"
echo "===================================================================================================="

"${PYTHON}" -u \
  src/imputation/common/compare_advisor_genes_across_experiments.py \
  --model "${MODEL}" \
  --model-label "${MODEL_LABEL}" \
  --combined-dir "${COMBINED_DIR}"

echo
echo "===================================================================================================="
echo "SUCCESS"
echo "Model: ${MODEL_LABEL}"
echo "Completed: $(date)"
echo "===================================================================================================="
