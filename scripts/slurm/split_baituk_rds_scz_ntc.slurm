#!/bin/bash
#SBATCH --job-name=baituk_rds_split
#SBATCH --output=logs/baituk_rds_split_%j.out
#SBATCH --error=logs/baituk_rds_split_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

R_ENV="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r"
RSCRIPT="${R_ENV}/bin/Rscript"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${PROJECT_ROOT}"
mkdir -p logs

export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "===================================================================================================="
echo "Batiuk RDS to separate SCZ and NTC H5AD"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "Memory: 200G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "R environment: ${R_ENV}"
echo "===================================================================================================="

echo
echo "STEP 1: Export sample-level H5AD files from the raw RDS"
echo

srun --cpu-bind=cores \
  "${RSCRIPT}" \
  src/snrnaseq_branch/export_baituk_rds_samples_to_h5ad.R

echo
echo "STEP 2: Concatenate sample H5AD files into SCZ and NTC outputs"
echo

srun --cpu-bind=cores \
  "${PYTHON}" -u \
  src/snrnaseq_branch/concat_baituk_samples_scz_ntc.py

echo
echo "===================================================================================================="
echo "Finished successfully: $(date)"
echo "===================================================================================================="
