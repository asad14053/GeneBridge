#!/bin/bash
#SBATCH --job-name=8667_genevar
#SBATCH --output=logs/gene_variation/%j.out
#SBATCH --error=logs/gene_variation/%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

export MPLBACKEND=Agg
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16

"${PYTHON}" \
    src/analysis/19_Br8667_OOF300_vs_34K_gene_variation.py

echo
echo "===== OUTPUTS ====="

ls -lh \
results/analysis/gene_variation/Br8667/
