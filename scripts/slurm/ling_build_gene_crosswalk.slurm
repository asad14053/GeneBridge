#!/bin/bash
#SBATCH --job-name=ling_xwalk
#SBATCH --output=logs/ling_harmonization/%j.out
#SBATCH --error=logs/ling_harmonization/%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

"${PYTHON}" \
    src/analysis/ling_harmonization/01_build_ling_h5_gene_crosswalk.py \
    --h5-root "${ROOT}/data/raw/nemo_mccarroll/counts" \
    --cache-dir "${ROOT}/data/reference/ling_gene_annotation_cache" \
    --out-dir "${ROOT}/results/analysis/ling_gene_harmonization"

echo
echo "===== OUTPUT FILES ====="

ls -lh \
    results/analysis/ling_gene_harmonization/
