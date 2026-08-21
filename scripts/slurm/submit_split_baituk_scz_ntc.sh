#!/bin/bash

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
SLURM_SCRIPT="${PROJECT_ROOT}/scripts/slurm/split_baituk_scz_ntc.slurm"
JOB_ID_FILE="${PROJECT_ROOT}/logs/latest_baituk_scz_ntc_job_id.txt"

cd "${PROJECT_ROOT}"
mkdir -p logs

# Always classify exact MB8 as SCZ
sed -i -E \
  's/--mb8-policy[[:space:]]+(error|ntc|scz)/--mb8-policy scz/' \
  "${SLURM_SCRIPT}"

bash -n "${SLURM_SCRIPT}"

JOB_ID=$(sbatch --parsable "${SLURM_SCRIPT}")

echo "${JOB_ID}" > "${JOB_ID_FILE}"

echo "Submitted Batiuk SCZ/NTC split job: ${JOB_ID}"

squeue -j "${JOB_ID}" \
  -o "%.18i %.30j %.10T %.12M %.12l %R"

echo
echo "Monitor:"
echo "tail -n 50 -F logs/baituk_scz_ntc_${JOB_ID}.out logs/baituk_scz_ntc_${JOB_ID}.err"
