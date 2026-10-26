#!/bin/bash

#SBATCH --job-name=xenium_N24_raw
#SBATCH --output=logs/xenium_N24_raw/%A_%a.out
#SBATCH --error=logs/xenium_N24_raw/%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --array=1-24%8

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

SCRIPT="src/xenium_branch/01_build_datastructure/01b_build_xenium_donor_raw.py"


DONORS=(
    Br1113
    Br1139
    Br2039
    Br2421
    Br2719
    Br5314
    Br5373
    Br5400
    Br5436
    Br5588
    Br5590
    Br5622
    Br5639
    Br5746
    Br5931
    Br5973
    Br6032
    Br6389
    Br6432
    Br6437
    Br6496
    Br8433
    Br8667
    Br8772
)


DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID - 1))]}"


cd "${PROJECT_ROOT}"

mkdir -p logs/xenium_N24_raw


export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"


echo "===================================================================================================="
echo "XENIUM DONOR RAW BUILD"
echo "===================================================================================================="
echo "Job ID       : ${SLURM_JOB_ID}"
echo "Array task   : ${SLURM_ARRAY_TASK_ID}"
echo "Donor        : ${DONOR}"
echo "Node         : $(hostname)"
echo "CPUs         : ${SLURM_CPUS_PER_TASK}"
echo "Started      : $(date)"
echo "===================================================================================================="


"${PYTHON}" -u "${SCRIPT}" \
    --sample "${DONOR}"


echo
echo "===================================================================================================="
echo "PASS"
echo "Donor     : ${DONOR}"
echo "Completed : $(date)"
echo "===================================================================================================="
