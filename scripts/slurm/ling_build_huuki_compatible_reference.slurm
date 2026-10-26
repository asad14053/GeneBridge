#!/bin/bash
#SBATCH --job-name=ling_ref
#SBATCH --output=logs/ling_reference_build/ling_ref_%j.out
#SBATCH --error=logs/ling_reference_build/ling_ref_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

# ============================================================
# FINAL LING REFERENCE
# ============================================================

REFERENCE_OUT="${ROOT}/data/processed/snrnaseq/ling_compatible/ling_snrna_reference_full_allgenes.h5ad"

cd "${ROOT}"

mkdir -p \
    logs/ling_reference_build \
    results/analysis/ling_reference_build \
    data/processed/snrnaseq/ling_compatible

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

WORK="/scratch/${USER}/ling_reference_${SLURM_JOB_ID}"

mkdir -p "${WORK}"

echo "============================================================"
echo "LING snRNA REFERENCE BUILD"
echo "============================================================"
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "Memory:      200G"
echo "Work dir:    ${WORK}"
echo "Final H5AD:  ${REFERENCE_OUT}"
echo "Started:     $(date)"
echo "============================================================"
echo


# ============================================================
# BUILD
# ============================================================

"${PYTHON}" \
    scripts/ling_harmonization/05_build_ling_huuki_compatible_reference.py \
    --root "${ROOT}" \
    --n-cells 77604 \
    --gene-prevalence 0.95 \
    --expected-genes 31491 \
    --seed 20261023 \
    --work-dir "${WORK}"


# ============================================================
# VERIFY FINAL REFERENCE
# ============================================================

echo
echo "============================================================"
echo "VERIFY FINAL OUTPUT"
echo "============================================================"

if [[ ! -f "${REFERENCE_OUT}" ]]; then
    echo "ERROR: Expected Ling reference was not created:"
    echo "${REFERENCE_OUT}"
    exit 1
fi

ls -lh "${REFERENCE_OUT}"


"${PYTHON}" - <<PY
import anndata as ad

p = "${REFERENCE_OUT}"

a = ad.read_h5ad(
    p,
    backed="r"
)

print()
print("Final Ling reference:")
print("Path:", p)
print("Shape:", a.shape)
print("OBS:", a.obs.shape)
print("VAR:", a.var.shape)
print("Layers:", list(a.layers.keys()))
print("OBSM:", list(a.obsm.keys()))
print("X dtype:", a.X.dtype)
print("var_names unique:", a.var_names.is_unique)

assert a.n_obs == 77604, a.shape
assert a.n_vars == 31491, a.shape
assert a.obs.shape[1] == 34, a.obs.shape
assert "logcounts" in a.layers
assert len(a.obsm.keys()) == 0
assert a.var_names.is_unique

a.file.close()

print()
print("REFERENCE VALIDATION: PASS")
PY


echo
echo "============================================================"
echo "BUILD COMPLETE"
echo "============================================================"
echo "Reference:"
echo "${REFERENCE_OUT}"
echo
echo "Finished: $(date)"
