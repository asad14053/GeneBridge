#!/bin/bash

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

cd "${ROOT}"

SLURM_SCRIPT="scripts/slurm/imputation/run_envi_N24_3fold_array.slurm"

LOG_DIR="${ROOT}/logs/envi_N24"

mkdir -p "${LOG_DIR}"


# ============================================================
# REFERENCES
# ============================================================

FULL_REF="${ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad"

BR6432_EX5="${ROOT}/data/processed/imputation_full/snrnaseq_target_references/Br6432/seq_data_huuki_snrna_Br6432.h5ad"

BR6432_EX53="${ROOT}/data/processed/imputation_full/snrnaseq_target_references/Br6432/seq_data_huuki_snrna_9brains_excluding_Br6432.h5ad"

BR8667_EX5="${ROOT}/data/processed/imputation_full/snrnaseq_target_references/Br8667/seq_data_huuki_snrna_Br8667.h5ad"

BR8667_EX53="${ROOT}/data/processed/imputation_full/snrnaseq_target_references/Br8667/seq_data_huuki_snrna_9brains_excluding_Br8667.h5ad"


# ============================================================
# VALIDATE REFERENCES
# ============================================================

for REF in \
    "${FULL_REF}" \
    "${BR6432_EX5}" \
    "${BR6432_EX53}" \
    "${BR8667_EX5}" \
    "${BR8667_EX53}"
do

    if [[ ! -f "${REF}" ]]; then
        echo "ERROR: missing reference:"
        echo "${REF}"
        exit 1
    fi

done


# ============================================================
# N24 TARGETS
# ============================================================

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


# ============================================================
# SUBMISSION MANIFEST
# ============================================================

MANIFEST="${LOG_DIR}/submission_manifest.tsv"

printf \
"job_id\tdonor\texperiment\treference\tfold_runs\n" \
> "${MANIFEST}"


echo "===================================================================================================="
echo "ENVI N24 SUBMISSION"
echo "===================================================================================================="

echo
echo "Primary:"
echo "  Ex-5.1 = all 24 donors × 3 folds"

echo
echo "Matched-reference sensitivity:"
echo "  Ex-5   = Br6432 + Br8667"
echo "  Ex-5.3 = Br6432 + Br8667"

echo
echo "Expected unique ENVI runs = 84"
echo


# ============================================================
# EX-5.1 — ALL 24 TARGETS
# ============================================================

echo "===================================================================================================="
echo "SUBMIT EX-5.1"
echo "===================================================================================================="

for DONOR in "${DONORS[@]}"; do

    JOBID=$(
        sbatch \
        --parsable \
        --job-name="envi_${DONOR}_ex51" \
        "${SLURM_SCRIPT}" \
        "${DONOR}" \
        "ex5_1" \
        "Experiment 5.1" \
        "${FULL_REF}"
    )

    echo \
    "Ex-5.1 | ${DONOR} | Job ${JOBID}"

    printf \
    "%s\t%s\t%s\t%s\t3\n" \
    "${JOBID}" \
    "${DONOR}" \
    "ex5_1" \
    "${FULL_REF}" \
    >> "${MANIFEST}"

done


# ============================================================
# EX-5 — MATCHED TARGET / MATCHED snRNA
# ============================================================

echo
echo "===================================================================================================="
echo "SUBMIT EX-5"
echo "===================================================================================================="

JOBID=$(
    sbatch \
    --parsable \
    --job-name="envi_Br6432_ex5" \
    "${SLURM_SCRIPT}" \
    "Br6432" \
    "ex5" \
    "Experiment 5" \
    "${BR6432_EX5}"
)

echo \
"Ex-5 | Br6432 | Job ${JOBID}"

printf \
"%s\t%s\t%s\t%s\t3\n" \
"${JOBID}" \
"Br6432" \
"ex5" \
"${BR6432_EX5}" \
>> "${MANIFEST}"


JOBID=$(
    sbatch \
    --parsable \
    --job-name="envi_Br8667_ex5" \
    "${SLURM_SCRIPT}" \
    "Br8667" \
    "ex5" \
    "Experiment 5" \
    "${BR8667_EX5}"
)

echo \
"Ex-5 | Br8667 | Job ${JOBID}"

printf \
"%s\t%s\t%s\t%s\t3\n" \
"${JOBID}" \
"Br8667" \
"ex5" \
"${BR8667_EX5}" \
>> "${MANIFEST}"


# ============================================================
# EX-5.3 — TARGET-EXCLUDED 9-BRAIN REFERENCE
# ============================================================

echo
echo "===================================================================================================="
echo "SUBMIT EX-5.3"
echo "===================================================================================================="

JOBID=$(
    sbatch \
    --parsable \
    --job-name="envi_Br6432_ex53" \
    "${SLURM_SCRIPT}" \
    "Br6432" \
    "ex5_3" \
    "Experiment 5.3" \
    "${BR6432_EX53}"
)

echo \
"Ex-5.3 | Br6432 | Job ${JOBID}"

printf \
"%s\t%s\t%s\t%s\t3\n" \
"${JOBID}" \
"Br6432" \
"ex5_3" \
"${BR6432_EX53}" \
>> "${MANIFEST}"


JOBID=$(
    sbatch \
    --parsable \
    --job-name="envi_Br8667_ex53" \
    "${SLURM_SCRIPT}" \
    "Br8667" \
    "ex5_3" \
    "Experiment 5.3" \
    "${BR8667_EX53}"
)

echo \
"Ex-5.3 | Br8667 | Job ${JOBID}"

printf \
"%s\t%s\t%s\t%s\t3\n" \
"${JOBID}" \
"Br8667" \
"ex5_3" \
"${BR8667_EX53}" \
>> "${MANIFEST}"


echo
echo "===================================================================================================="
echo "SUBMISSION COMPLETE"
echo "===================================================================================================="

echo "Job arrays submitted : 28"
echo "Fold runs submitted  : 84"
echo
echo "Manifest:"
echo "${MANIFEST}"
