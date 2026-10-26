#!/bin/bash
#SBATCH --job-name=enviPI12
#SBATCH --output=logs/envi_PI_12NTC/%A_%a.out
#SBATCH --error=logs/envi_PI_12NTC/%A_%a.err
#SBATCH --array=0-11%12
#SBATCH --time=06:00:00
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

# Avoid hidden nested threading.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "============================================================"
echo "ENVI PI FIGURE"
echo "BR=${BR}"
echo "TASK=${SLURM_ARRAY_TASK_ID}"
echo "============================================================"

"${PYTHON}" \
    src/analysis/08_envi_fulltranscriptome_PI_figure.py

echo
echo "${BR} PASS"
