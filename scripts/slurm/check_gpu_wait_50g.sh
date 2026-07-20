#!/bin/bash

set -euo pipefail

ACCOUNT="${1:-ai-gpu}"
MEM="${2:-50G}"
TIME_LIMIT="${3:-12:00:00}"
CPUS="${4:-8}"

PARTITIONS=(
  "b200-8-gm1432-c192-m2048"
  "rp6b-1-gm96-c8-m64"
  "rp6b-1-gm96-c8-m64-bk"
  "rp6b-8-gm768-c192-m2048"
  "l4-4-gm96-c48-m192"
  "l4-4-gm96-c48-m192-bk"
)

echo "============================================================"
echo "GPU wait-time check"
echo "Account:     ${ACCOUNT}"
echo "Memory:      ${MEM}"
echo "CPUs:        ${CPUS}"
echo "GPUs:        1"
echo "Time limit:  ${TIME_LIMIT}"
echo "============================================================"
echo

for PARTITION in "${PARTITIONS[@]}"; do
    echo "------------------------------------------------------------"
    echo "Partition: ${PARTITION}"
    echo "------------------------------------------------------------"

    echo "[sinfo]"
    sinfo -h -p "${PARTITION}" -o "partition=%P state=%t nodes=%D gres=%G mem_mb=%m cpus=%C" 2>/dev/null || true

    echo
    echo "[current queue summary]"
    squeue -h -p "${PARTITION}" -o "%T" 2>/dev/null | sort | uniq -c || true

    TMP_SCRIPT=$(mktemp)

    cat > "${TMP_SCRIPT}" <<EOT
#!/bin/bash
#SBATCH --job-name=gpu_wait_test
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --gpus=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

sleep 60
EOT

    echo
    echo "[sbatch --test-only estimate]"
    sbatch --test-only "${TMP_SCRIPT}" 2>&1 | sed 's/^/  /' || true

    rm -f "${TMP_SCRIPT}"

    echo
done

echo "============================================================"
echo "Done."
echo "Pick the partition with the earliest estimated start time."
echo "If a partition says inactive/drain/unavailable, skip it."
echo "============================================================"
