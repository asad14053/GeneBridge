#!/bin/bash
#SBATCH --job-name=8667_morancmp
#SBATCH --output=logs/advisor_moran_cosine/%j.out
#SBATCH --error=logs/advisor_moran_cosine/%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

export MPLBACKEND=Agg
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

"${PYTHON}" \
    src/analysis/21_Br8667_advisor_moran_cosine_comparison.py

echo
echo "===== OUTPUT CHECK ====="

ls -lh \
results/analysis/three_source_validation/Br8667/advisor_moran_cosine_validation/
