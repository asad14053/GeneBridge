#!/bin/bash
#SBATCH --job-name=gene_harmonize
#SBATCH --output=logs/ling_harmonization/%j.out
#SBATCH --error=logs/ling_harmonization/%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

"${PYTHON}" \
scripts/ling_harmonization/02_harmonize_ling_huuki_kwon_genes.py \
--ling-h5-root "${ROOT}/data/raw/nemo_mccarroll/counts" \
--huuki-h5ad "${ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad" \
--kwon-h5ad "${ROOT}/data/processed/imputation_full/ex1_ntc/targets/spatial_data_xenium_Br8667_ex1_ntc.h5ad" \
--cache-dir "${ROOT}/data/reference/ling_gene_annotation_cache" \
--out-dir "${ROOT}/results/analysis/ling_huuki_kwon_gene_harmonization" \
--expected-ling-files 87

echo
echo "===== OUTPUT ====="

ls -lh \
"${ROOT}/results/analysis/ling_huuki_kwon_gene_harmonization"
