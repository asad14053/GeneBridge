#!/bin/bash

#SBATCH --job-name=merge_xenium_N24
#SBATCH --output=logs/merge_xenium_N24_%j.out
#SBATCH --error=logs/merge_xenium_N24_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G

set -euo pipefail

cd /beegfs/labs/hulab/projects/mjabin/GeneBridge

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"${PYTHON}" -u \
src/xenium_branch/01_build_datastructure/02_merge_xenium_N24_raw.py
