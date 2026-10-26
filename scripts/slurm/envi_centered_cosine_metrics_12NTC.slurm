#!/bin/bash
#SBATCH --job-name=enviCentCos
#SBATCH --output=logs/envi_centered_cosine_12NTC/metrics_%A_%a.out
#SBATCH --error=logs/envi_centered_cosine_12NTC/metrics_%A_%a.err
#SBATCH --array=0-11%12
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

DONORS=(
    Br1113
    Br2719
    Br5314
    Br5400
    Br5436
    Br5622
    Br5639
    Br5931
    Br6389
    Br6432
    Br8433
    Br8667
)

export BR="${DONORS[$SLURM_ARRAY_TASK_ID]}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "============================================================"
echo "CENTERED COSINE METRICS"
echo "BR=${BR}"
echo "TASK=${SLURM_ARRAY_TASK_ID}"
echo "============================================================"

"${PYTHON}" \
    src/analysis/06_envi_fulltranscriptome_crosstech.py

echo
echo "${BR} METRICS PASS"
