#!/bin/bash
#SBATCH --job-name=N24_gene_folds
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=1-24%12
#SBATCH --output=logs/gene_folds_N24/%A_%a.out
#SBATCH --error=logs/gene_folds_N24/%A_%a.err

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

INPUT="${ROOT}/data/processed/xenium/xenium_N24_layer_celltype_annotated.h5ad"

CANONICAL="${ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100/gene_metrics_and_fold_assignment.h5ad"

OUT="${ROOT}/data/processed/imputation_full/gene_folds_200_100"

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

DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "============================================================"
echo "SLURM task : ${SLURM_ARRAY_TASK_ID}"
echo "Donor      : ${DONOR}"
echo "Input      : ${INPUT}"
echo "Canonical  : ${CANONICAL}"
echo "Output     : ${OUT}/${DONOR}"
echo "============================================================"

mkdir -p "${OUT}/${DONOR}"

"${PYTHON}" \
    "${ROOT}/src/imputation/apply_canonical_gene_folds_N24.py" \
    --input "${INPUT}" \
    --canonical-folds "${CANONICAL}" \
    --donor "${DONOR}" \
    --output-root "${OUT}"

echo
echo "DONE: ${DONOR}"
