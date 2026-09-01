#!/bin/bash
#SBATCH --job-name=xenium_download
#SBATCH --output=logs/download_xenium/%A_%a.out
#SBATCH --error=logs/download_xenium/%A_%a.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --array=1-24%24

set -euo pipefail

cd /beegfs/labs/hulab/projects/mjabin/GeneBridge

PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

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
echo "XENIUM PARALLEL DOWNLOAD"
echo "============================================================"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Donor      : ${DONOR}"
echo "Node       : $(hostname)"
echo "CPUs       : ${SLURM_CPUS_PER_TASK}"
echo "Memory     : 200G"
echo "Started    : $(date)"
echo "============================================================"

${PYTHON} -u \
src/data/download_24_common_xenium_visium.py \
--donor "${DONOR}" \
--include-heavy

echo
echo "============================================================"
echo "FINISHED"
echo "Donor   : ${DONOR}"
echo "Finished: $(date)"
echo "============================================================"
