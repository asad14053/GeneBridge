#!/bin/bash
#SBATCH --job-name=linglib
#SBATCH --output=logs/ling_reference_build/linglib_%A_%a.out
#SBATCH --error=logs/ling_reference_build/linglib_%A_%a.err
#SBATCH --array=0-86%24
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"${PYTHON}" \
scripts/ling_harmonization/07_process_ling_library.py \
--task-id "${SLURM_ARRAY_TASK_ID}"
