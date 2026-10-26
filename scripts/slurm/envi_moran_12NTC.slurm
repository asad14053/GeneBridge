#!/bin/bash
#SBATCH --job-name=enviMoran
#SBATCH --output=logs/envi_moran_12NTC/%A_%a.out
#SBATCH --error=logs/envi_moran_12NTC/%A_%a.err
#SBATCH --array=0-11%12
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4

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

BR="${DONORS[$SLURM_ARRAY_TASK_ID]}"

export BR
export MORAN_K=6

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

echo "============================================================"
echo "BR=${BR}"
echo "MORAN_K=${MORAN_K}"
echo "============================================================"

"${PYTHON}" \
    src/analysis/09_envi_moran_visium_vs_imputed.py

echo
echo "${BR} PASS"
