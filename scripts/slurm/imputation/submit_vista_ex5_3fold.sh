
#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
cd "${PROJECT_ROOT}"

mkdir -p logs

ARRAY_JOB_ID=$(
    sbatch --parsable \
        scripts/slurm/imputation/run_vista_ex5_3fold_array.slurm
)

AGGREGATE_JOB_ID=$(
    sbatch --parsable \
        --dependency="afterok:${ARRAY_JOB_ID}" \
        scripts/slurm/imputation/aggregate_vista_ex5.slurm
)

echo "Submitted VISTA 3-fold array: ${ARRAY_JOB_ID}"
echo "Submitted dependent aggregation: ${AGGREGATE_JOB_ID}"
echo
echo "Monitor:"
echo "  watch -n 10 \"squeue -j ${ARRAY_JOB_ID},${AGGREGATE_JOB_ID} -o '%.12i %.30j %.10T %.12M %.12l %R'\""
