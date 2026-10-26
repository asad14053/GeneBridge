#!/bin/bash
#SBATCH --job-name=8667_genewise
#SBATCH --output=logs/advisor_genewise/%j.out
#SBATCH --error=logs/advisor_genewise/%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
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
    src/analysis/20_Br8667_advisor_genewise_OOF_figure.py

echo
echo "===== OUTPUT CHECK ====="

ls -lh \
results/analysis/three_source_validation/Br8667/advisor_genewise_validation/
