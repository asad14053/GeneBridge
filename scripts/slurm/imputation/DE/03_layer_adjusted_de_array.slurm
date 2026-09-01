#!/bin/bash
#SBATCH --job-name=layer_adj_de
#SBATCH --output=logs/imputation/DE/layer_adjusted/%A_%a.out
#SBATCH --error=logs/imputation/DE/layer_adjusted/%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --array=1-2%2

set -euo pipefail

cd /beegfs/labs/hulab/projects/mjabin/GeneBridge

RSCRIPT=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript

case "${SLURM_ARRAY_TASK_ID}" in

    1)
        DATASET="original"
        ;;

    2)
        DATASET="envi_measured"
        ;;

    *)
        echo "Invalid array task."
        exit 1
        ;;

esac


echo "=========================================="
echo "Dataset     : ${DATASET}"
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Array task  : ${SLURM_ARRAY_TASK_ID}"
echo "Node        : $(hostname)"
echo "CPUs        : ${SLURM_CPUS_PER_TASK}"
echo "Memory      : 200G"
echo "=========================================="


${RSCRIPT} \
src/imputation/DE/paper_concordance/03_layer_adjusted_de.R \
"${DATASET}"
