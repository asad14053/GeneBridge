#!/bin/bash
#SBATCH --job-name=oof_gene_plot
#SBATCH --output=logs/oof_gene_plot_%j.out
#SBATCH --error=logs/oof_gene_plot_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/common/plot_oof_advisor_genes.py"

EXPERIMENT="${1:?Missing experiment key}"
EXPERIMENT_LABEL="${2:?Missing experiment label}"
MODEL="${3:?Missing model key}"
MODEL_LABEL="${4:?Missing model label}"
COMBINED_DIR="${5:-auto}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

cd "${PROJECT_ROOT}"

mkdir -p logs

echo "===================================================================================================="
echo "Fixed advisor-gene OOF spatial comparison"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Memory: 160G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Experiment: ${EXPERIMENT}"
echo "Model: ${MODEL}"
echo "Combined directory: ${COMBINED_DIR}"
echo "===================================================================================================="

"${PYTHON}" -m py_compile \
  "${SCRIPT}" \
  "${PROJECT_ROOT}/src/imputation/common/benchmark_evaluation.py"

"${PYTHON}" -u "${SCRIPT}" \
  --experiment "${EXPERIMENT}" \
  --experiment-label "${EXPERIMENT_LABEL}" \
  --model "${MODEL}" \
  --model-label "${MODEL_LABEL}" \
  --combined-dir "${COMBINED_DIR}"

echo
echo "Completed: $(date)"
