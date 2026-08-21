#!/bin/bash
#SBATCH --job-name=gimvi_exp_analysis
#SBATCH --output=logs/gimvi_exp_analysis_%j.out
#SBATCH --error=logs/gimvi_exp_analysis_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

cd "${PROJECT_ROOT}"

mkdir -p logs

export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

OUTPUT_DIR="${PROJECT_ROOT}/outputs/imputation_beta/Br8667/comparisons/gimvi"

echo "===================================================================================================="
echo "gimVI experiment analysis"
echo "===================================================================================================="
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Output: ${OUTPUT_DIR}"
echo

echo "===================================================================================================="
echo "STEPS 1-3: paired t-test + TOST + uncertainty plots"
echo "===================================================================================================="

"${PYTHON}" -u \
  src/imputation/common/paired_tests_imputation_experiments.py \
  --project-root "${PROJECT_ROOT}" \
  --model gimvi \
  --model-label "gimVI" \
  --aggregate-subdir combined_v2 \
  --output-dir "${OUTPUT_DIR}/paired_experiment_tests" \
  --seed 8667 \
  --n-draws 10000

echo
echo "===================================================================================================="
echo "STEP 4: advisor 11-gene panel"
echo "===================================================================================================="

"${PYTHON}" -u \
  src/imputation/common/advisor_gene_panel_imputation.py \
  --project-root "${PROJECT_ROOT}" \
  --model gimvi \
  --model-label "gimVI" \
  --aggregate-subdir combined_v2 \
  --output-dir "${OUTPUT_DIR}/advisor_11_genes"

echo
echo "===================================================================================================="
echo "SUCCESS"
echo "Completed: $(date)"
echo "===================================================================================================="
