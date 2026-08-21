#!/bin/bash
#SBATCH --job-name=baituk_scz_ntc
#SBATCH --output=logs/baituk_scz_ntc_%j.out
#SBATCH --error=logs/baituk_scz_ntc_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

INPUT="${PROJECT_ROOT}/data/processed/snrnaseq/Baituk/baituk_snrna_reference_full_allgenes.h5ad"

SCZ_OUTPUT="${PROJECT_ROOT}/data/processed/snrnaseq/Baituk/baituk_snrna_reference_SCZ_full_allgenes.h5ad"

NTC_OUTPUT="${PROJECT_ROOT}/data/processed/snrnaseq/Baituk/baituk_snrna_reference_NTC_full_allgenes.h5ad"

CHUNK_DIR="${PROJECT_ROOT}/data/processed/snrnaseq/Baituk/baituk_scz_ntc_chunks"

cd "${PROJECT_ROOT}"
mkdir -p logs

export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

echo "============================================================"
echo "Chunked Batiuk SCZ/NTC split"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "Input: ${INPUT}"
echo "Chunk directory: ${CHUNK_DIR}"
echo "============================================================"

srun --cpu-bind=cores \
  "${PYTHON}" -u \
  src/snrnaseq_branch/split_baituk_scz_ntc_chunked.py \
  --input "${INPUT}" \
  --scz-output "${SCZ_OUTPUT}" \
  --ntc-output "${NTC_OUTPUT}" \
  --chunk-dir "${CHUNK_DIR}" \
  --chunk-size 5000 \
  --retries 3

echo "Finished: $(date)"
