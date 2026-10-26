#!/bin/bash
#SBATCH --job-name=t23_pb
#SBATCH --array=1-6%6
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/task2_3_pb/%A_%a.out
#SBATCH --error=logs/task2_3_pb/%A_%a.err

set -euo pipefail

cd /beegfs/labs/hulab/projects/mjabin/GeneBridge

MODELS=(VISTA gimVI ENVI TransImp Tangram SpaGE)
MODEL="${MODELS[$((SLURM_ARRAY_TASK_ID-1))]}"

PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

echo "================================================================================"
echo "TASK 2.3 — OOF DONOR-DOMAIN PSEUDOBULK"
echo "================================================================================"
echo "Task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Model   : ${MODEL}"
echo "Python  : ${PYTHON}"
echo "Started : $(date)"
echo

"${PYTHON}" \
src/imputation/DE/task2_3_ntc_dysregulation/01_build_model_oof_layer_pseudobulk.py \
--model "${MODEL}"

echo
echo "================================================================================"
echo "MODEL: ${MODEL}"
echo "SLURM STATUS: PASS"
echo "Finished: $(date)"
echo "================================================================================"
