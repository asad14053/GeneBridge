#!/bin/bash
#SBATCH --job-name=enviXT8667
#SBATCH --output=logs/envi_full_crosstech/Br8667_%j.out
#SBATCH --error=logs/envi_full_crosstech/Br8667_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

export BR="Br8667"

# We parallelize explicitly across gene chunks.
# Prevent nested BLAS/OpenMP oversubscription.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export N_JOBS="${SLURM_CPUS_PER_TASK}"

"${PYTHON}" \
    src/analysis/06_envi_fulltranscriptome_crosstech.py
