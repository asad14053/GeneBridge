#!/bin/bash
#SBATCH --job-name=lingmerge
#SBATCH --output=logs/ling_reference_build/lingmerge_%j.out
#SBATCH --error=logs/ling_reference_build/lingmerge_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

WORK="/scratch/${USER}/ling_reference_merge_${SLURM_JOB_ID}"

mkdir -p "${WORK}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"${PYTHON}" \
scripts/ling_harmonization/08_merge_ling_reference.py \
--work-dir "${WORK}"
