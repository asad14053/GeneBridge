#!/bin/bash
#SBATCH --job-name=envi_layer_pb
#SBATCH --output=logs/imputation/DE/envi_layer_pseudobulk/%A_%a.out
#SBATCH --error=logs/imputation/DE/envi_layer_pseudobulk/%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --array=1-23%23

set -euo pipefail

cd /beegfs/labs/hulab/projects/mjabin/GeneBridge

PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "=========================================="
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Array task  : ${SLURM_ARRAY_TASK_ID}"
echo "Node        : $(hostname)"
echo "CPUs        : ${SLURM_CPUS_PER_TASK}"
echo "Memory      : 200G"
echo "=========================================="

${PYTHON} \
src/imputation/DE/paper_concordance/02a_build_envi_layer_pseudobulk_donor.py \
--index "${SLURM_ARRAY_TASK_ID}"
