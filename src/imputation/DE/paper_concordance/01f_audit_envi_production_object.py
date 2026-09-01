#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import anndata as ad
import scipy.sparse as sp


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

DONOR = "Br1113"

TARGET = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/targets"
    / f"spatial_data_xenium_{DONOR}_ex1_ntc.h5ad"
)

ENVI = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/envi"
    / DONOR
    / f"spatial_data_xenium_{DONOR}_ENVI_full_transcriptome.h5ad"
)


def section(x):
    print("\n" + "=" * 110)
    print(x)
    print("=" * 110)


def dense(x):
    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


section("FILES")

print("Target:", TARGET)
print("Target exists:", TARGET.exists())

print("ENVI:", ENVI)
print("ENVI exists:", ENVI.exists())


section("LOAD BACKED OBJECTS")

target = ad.read_h5ad(
    TARGET,
    backed="r",
)

envi = ad.read_h5ad(
    ENVI,
    backed="r",
)

print("Target shape:", target.shape)
print("ENVI shape  :", envi.shape)


section("TARGET STRUCTURE")

print("layers:")
print(list(target.layers.keys()))

print("\nobs columns:")
print(target.obs.columns.tolist())

print("\nvar columns:")
print(target.var.columns.tolist())


section("ENVI STRUCTURE")

print("layers:")
print(list(envi.layers.keys()))

print("\nobs columns:")
print(envi.obs.columns.tolist())

print("\nvar columns:")
print(envi.var.columns.tolist())


section("SEARCH ENVI METADATA FOR MEASURED / IMPUTED / OOF INFORMATION")

interesting_cols = [
    c for c in envi.var.columns
    if any(
        k in c.lower()
        for k in [
            "source",
            "measured",
            "imputed",
            "oof",
            "predict",
            "panel",
        ]
    )
]

print("Interesting var columns:")
print(interesting_cols)

for c in interesting_cols:

    print(f"\nCOLUMN: {c}")

    try:
        print(
            envi.var[c]
            .astype(str)
            .value_counts(dropna=False)
            .head(30)
            .to_string()
        )
    except Exception as e:
        print("ERROR:", repr(e))


section("SEARCH LAYER NAMES FOR OOF / PREDICTION")

possible_prediction_layers = [
    x for x in envi.layers.keys()
    if any(
        k in x.lower()
        for k in [
            "oof",
            "predict",
            "imput",
            "heldout",
        ]
    )
]

print(
    "Prediction-like ENVI layers:",
    possible_prediction_layers
)


section("GENE UNIVERSE")

target_genes = set(
    target.var_names.astype(str)
)

envi_genes = set(
    envi.var_names.astype(str)
)

shared = sorted(
    target_genes & envi_genes
)

print("Target genes:", len(target_genes))
print("ENVI genes  :", len(envi_genes))
print("Shared      :", len(shared))

print(
    "All target 300 genes present in ENVI:",
    len(shared) == target.n_vars
)


section("CELL UNIVERSE")

target_cells = set(
    target.obs_names.astype(str)
)

envi_cells = set(
    envi.obs_names.astype(str)
)

shared_cells = sorted(
    target_cells & envi_cells
)

print("Target cells:", len(target_cells))
print("ENVI cells  :", len(envi_cells))
print("Shared cells:", len(shared_cells))

print(
    "Same cell universe:",
    target_cells == envi_cells
)


section("DIRECT TEST: ORIGINAL 300 vs ENVI SAME 300")

# Use a manageable sample first.
test_cells = shared_cells[:100]
test_genes = shared[:100]

target_cell_idx = target.obs_names.get_indexer(
    test_cells
)

envi_cell_idx = envi.obs_names.get_indexer(
    test_cells
)

target_gene_idx = target.var_names.get_indexer(
    test_genes
)

envi_gene_idx = envi.var_names.get_indexer(
    test_genes
)


# Prefer raw count layer in original Xenium.
if "counts" in target.layers.keys():

    original = target.layers["counts"][
        target_cell_idx,
        target_gene_idx
    ]

    original_source = "target.layers['counts']"

else:

    original = target.X[
        target_cell_idx,
        target_gene_idx
    ]

    original_source = "target.X"


envi_values = envi.X[
    envi_cell_idx,
    envi_gene_idx
]


original = dense(original).astype(float)
envi_values = dense(envi_values).astype(float)


print("Original source:", original_source)
print("ENVI source    : ENVI.X")

print("Compared values:", original.size)

print(
    "Exact values:",
    np.sum(
        original == envi_values
    ),
)

print(
    "Exact fraction:",
    np.mean(
        original == envi_values
    ),
)

print(
    "Maximum absolute difference:",
    np.max(
        np.abs(
            original - envi_values
        )
    ),
)


section("IMPORTANT INTERPRETATION")

if np.array_equal(
    original,
    envi_values,
):

    print(
        """
The tested ENVI values for the original Xenium genes are exactly
the same as the measured Xenium values.

Therefore these genes are being PRESERVED in the production
full-transcriptome ENVI object.

They are not evidence of OOF prediction.
"""
    )

else:

    print(
        """
The tested original genes differ between target and ENVI.X.

We must determine whether these values are predictions,
normalization-transformed values, or another representation.
"""
    )


if len(possible_prediction_layers) == 0:

    print(
        """
No obvious OOF / held-out / prediction layer was found in this
production ENVI object.
"""
    )

else:

    print(
        "\nPotential prediction layers:",
        possible_prediction_layers,
    )


target.file.close()
envi.file.close()

print("\nFINAL STATUS: AUDIT COMPLETE")
