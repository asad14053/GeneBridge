#!/usr/bin/env python3

from pathlib import Path
import re

import anndata as ad
import pandas as pd


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

META = (
    ROOT
    / "data/metadata/xenium_DE_metadata_23.csv"
)

ANNOTATED = (
    ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad"
)

TARGET_DIRS = [
    ROOT / "data/processed/imputation_full/ex1_ntc/targets",
    ROOT / "data/processed/imputation_full/ex2_scz/targets",
]

OUT = (
    ROOT
    / "outputs/imputation_full/DE/"
      "paper_concordance/layer_audit"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

EXPECTED_LAYERS = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WMtz",
    "WM",
]

EXPECTED_SPD_TO_LAYER = {
    "spd07": "L1/M",
    "spd06": "L2/3",
    "spd02": "L3/4",
    "spd05": "L5",
    "spd03": "L6",
    "spd01": "WMtz",
    "spd04": "WM",
}

MIN_NCELLS = 10


def donor_from_filename(path):
    m = re.search(r"(Br\d+)", path.name)

    if m is None:
        raise RuntimeError(
            f"Could not identify donor from {path.name}"
        )

    return m.group(1)


print("=" * 100)
print("PAPER-CONCORDANCE STEP 1")
print("23-DONOR × 7-SPATIAL-LAYER COVERAGE AUDIT")
print("=" * 100)


# =============================================================================
# Canonical 23-donor metadata
# =============================================================================

meta = pd.read_csv(META)

required = {
    "BrNum",
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "run_date",
}

missing = required - set(meta.columns)

if missing:
    raise RuntimeError(
        f"Missing canonical metadata columns: {sorted(missing)}"
    )

meta["BrNum"] = (
    meta["BrNum"]
    .astype(str)
    .str.strip()
)

assert len(meta) == 23
assert meta["BrNum"].nunique() == 23
assert "Br6432" not in set(meta["BrNum"])

expected_donors = set(meta["BrNum"])

print("\nCanonical cohort:")
print("  donors:", len(meta))
print("  NTC:", (meta["Dx"] == "NTC").sum())
print("  SCZ:", (meta["Dx"] == "SCZ").sum())
print("  Br6432 excluded:", "Br6432" not in expected_donors)


# =============================================================================
# Load annotated N24 object
# =============================================================================

print("\nLoading annotated N24 Xenium object...")

ann = ad.read_h5ad(
    ANNOTATED,
    backed="r",
)

required_obs = {
    "BrNum",
    "layer_annotation",
    "predictions_smooth",
    "cell_type_annotation",
}

missing_obs = required_obs - set(ann.obs.columns)

if missing_obs:
    raise RuntimeError(
        f"Annotated Xenium object missing: {sorted(missing_obs)}"
    )


annotation = ann.obs[
    [
        "BrNum",
        "layer_annotation",
        "predictions_smooth",
        "cell_type_annotation",
    ]
].copy()

annotation["BrNum"] = (
    annotation["BrNum"]
    .astype(str)
    .str.strip()
)

annotation["layer_annotation"] = (
    annotation["layer_annotation"]
    .astype(str)
    .str.strip()
)

annotation["predictions_smooth"] = (
    annotation["predictions_smooth"]
    .astype(str)
    .str.strip()
)


print("Annotated object shape:", ann.shape)
print(
    "Annotated donors:",
    annotation["BrNum"].nunique(),
)


# =============================================================================
# Validate SpD → layer mapping
# =============================================================================

mapping = (
    annotation[
        [
            "predictions_smooth",
            "layer_annotation",
        ]
    ]
    .drop_duplicates()
    .sort_values("predictions_smooth")
)


print("\nObserved SpD → layer mapping:")
print(mapping.to_string(index=False))


for spd, layer in EXPECTED_SPD_TO_LAYER.items():

    rows = mapping.loc[
        mapping["predictions_smooth"] == spd,
        "layer_annotation",
    ].tolist()

    if rows != [layer]:
        raise RuntimeError(
            f"Unexpected mapping for {spd}: "
            f"observed={rows}, expected={layer}"
        )


# =============================================================================
# Locate target H5ADs
# =============================================================================

files = []

for directory in TARGET_DIRS:

    if not directory.exists():
        raise FileNotFoundError(directory)

    files.extend(
        sorted(directory.glob("*.h5ad"))
    )


print("\nTarget H5AD files:", len(files))

if len(files) != 23:
    raise RuntimeError(
        f"Expected 23 target H5ADs, found {len(files)}"
    )


# =============================================================================
# Match targets to annotations and count donor × layer cells
# =============================================================================

coverage_rows = []
match_rows = []
exclusion_rows = []

observed_donors = set()


for i, path in enumerate(files, start=1):

    donor = donor_from_filename(path)

    print(
        f"\n[{i:02d}/23] "
        f"{donor} | {path.name}"
    )

    if donor not in expected_donors:
        raise RuntimeError(
            f"{donor} is not in canonical 23-donor metadata."
        )

    if donor in observed_donors:
        raise RuntimeError(
            f"Duplicate target donor: {donor}"
        )

    observed_donors.add(donor)


    target = ad.read_h5ad(
        path,
        backed="r",
    )


    target_ids = pd.Index(
        target.obs_names.astype(str)
    )


    ann_donor = annotation.loc[
        annotation["BrNum"] == donor
    ].copy()


    ann_ids = pd.Index(
        ann_donor.index.astype(str)
    )


    overlap = target_ids.intersection(
        ann_ids
    )


    n_target = len(target_ids)
    n_ann = len(ann_ids)
    n_overlap = len(overlap)


    target_fraction = (
        n_overlap / n_target
        if n_target > 0
        else 0
    )

    annotated_fraction = (
        n_overlap / n_ann
        if n_ann > 0
        else 0
    )


    print("  target cells:", n_target)
    print("  annotated cells:", n_ann)
    print("  matching IDs:", n_overlap)
    print(
        "  target overlap fraction:",
        f"{target_fraction:.6f}"
    )
    print(
        "  annotated overlap fraction:",
        f"{annotated_fraction:.6f}"
    )


    target_only = target_ids.difference(ann_ids)
    annotated_only = ann_ids.difference(target_ids)

    print("  target-only cells:", len(target_only))
    print("  annotated-only cells:", len(annotated_only))


    # This is the unsafe direction:
    # a cell has an annotation but no corresponding expression target.
    if len(annotated_only) > 0:
        raise RuntimeError(
            f"{donor}: {len(annotated_only)} annotated cells "
            "are missing from the expression target."
        )


    # Safe/documented exclusion:
    # expression exists, but no paper-compatible layer annotation.
    if len(target_only) > 0:

        print(
            f"  WARNING: excluding {len(target_only)} "
            "target cells without layer annotation."
        )

        for cell_id in target_only:
            exclusion_rows.append(
                {
                    "BrNum": donor,
                    "cell_id": str(cell_id),
                    "reason": "missing_layer_annotation",
                    "target_file": str(path),
                }
            )


    # Use only cells shared between expression target
    # and annotated N24 object.
    matched_ids = target_ids.intersection(ann_ids)

    ann_match = annotation.loc[
        matched_ids
    ].copy()


    if not (
        ann_match["BrNum"] == donor
    ).all():
        raise RuntimeError(
            f"{donor}: BrNum mismatch after cell-ID join."
        )


    layer_counts = (
        ann_match["layer_annotation"]
        .value_counts()
        .reindex(
            EXPECTED_LAYERS,
            fill_value=0,
        )
    )


    print("  layer counts:")
    print(layer_counts.to_string())


    for layer, n_cells in layer_counts.items():

        coverage_rows.append(
            {
                "BrNum": donor,
                "layer_annotation": layer,
                "n_cells": int(n_cells),
                "passes_min10":
                    int(n_cells) >= MIN_NCELLS,
            }
        )


    match_rows.append(
        {
            "BrNum": donor,
            "target_cells": n_target,
            "annotated_cells": n_ann,
            "matching_cells": n_overlap,
            "target_overlap_fraction":
                target_fraction,
            "annotated_overlap_fraction":
                annotated_fraction,
            "target_only_unannotated_cells":
                len(target_only),
            "annotated_only_missing_target_cells":
                len(annotated_only),
            "target_file": str(path),
        }
    )


    target.file.close()


ann.file.close()


# =============================================================================
# Cohort validation
# =============================================================================

if observed_donors != expected_donors:

    raise RuntimeError(
        "Target donor set differs from canonical metadata.\n"
        f"Missing: {sorted(expected_donors - observed_donors)}\n"
        f"Unexpected: {sorted(observed_donors - expected_donors)}"
    )


coverage = pd.DataFrame(
    coverage_rows
)

matches = pd.DataFrame(
    match_rows
)

exclusions = pd.DataFrame(
    exclusion_rows,
    columns=[
        "BrNum",
        "cell_id",
        "reason",
        "target_file",
    ],
)


# =============================================================================
# Validate layers
# =============================================================================

observed_layers = set(
    coverage.loc[
        coverage["n_cells"] > 0,
        "layer_annotation",
    ]
)

expected_layers = set(
    EXPECTED_LAYERS
)


print("\n" + "=" * 100)
print("SPATIAL LAYER VALIDATION")
print("=" * 100)

print(
    "Observed layers:",
    sorted(observed_layers),
)

print(
    "Expected layers:",
    EXPECTED_LAYERS,
)


if observed_layers != expected_layers:
    raise RuntimeError(
        "Observed spatial layers do not match "
        "expected paper-compatible layers."
    )


# =============================================================================
# Coverage matrices
# =============================================================================

wide = (
    coverage
    .pivot(
        index="BrNum",
        columns="layer_annotation",
        values="n_cells",
    )
    .reindex(
        index=sorted(expected_donors),
        columns=EXPECTED_LAYERS,
    )
)


summary = (
    coverage
    .groupby("layer_annotation")
    ["n_cells"]
    .agg(
        donors_present=lambda x:
            int((x > 0).sum()),

        donors_min10=lambda x:
            int((x >= MIN_NCELLS).sum()),

        min="min",
        median="median",
        max="max",
    )
    .reindex(EXPECTED_LAYERS)
)


bad = coverage.loc[
    ~coverage["passes_min10"]
].copy()


# =============================================================================
# Print results
# =============================================================================

print("\n" + "=" * 100)
print("23 DONOR × 7 LAYER CELL COUNTS")
print("=" * 100)

print(
    wide.to_string()
)


print("\n" + "=" * 100)
print("LAYER COVERAGE SUMMARY")
print("=" * 100)

print(
    summary.to_string()
)


print("\n" + "=" * 100)
print("DONOR × LAYER COMBINATIONS <10 CELLS")
print("=" * 100)

if bad.empty:
    print("NONE")
else:
    print(
        bad.to_string(index=False)
    )


# =============================================================================
# Save outputs
# =============================================================================

coverage.to_csv(
    OUT
    / "expected_23x7_layer_cell_counts_long.csv",
    index=False,
)

wide.to_csv(
    OUT
    / "expected_23x7_layer_cell_counts_matrix.csv"
)

summary.to_csv(
    OUT
    / "layer_coverage_summary.csv"
)

matches.to_csv(
    OUT
    / "target_annotation_cell_id_match_summary.csv",
    index=False,
)

mapping.to_csv(
    OUT
    / "spd_to_layer_mapping.csv",
    index=False,
)

exclusions.to_csv(
    OUT
    / "target_cells_excluded_missing_layer_annotation.csv",
    index=False,
)

bad.to_csv(
    OUT
    / "donor_layer_below_min10.csv",
    index=False,
)


# =============================================================================
# Final
# =============================================================================

n_total = len(coverage)

n_pass = int(
    coverage["passes_min10"].sum()
)

n_fail = (
    n_total
    - n_pass
)


print("\n" + "=" * 100)
print("FINAL AUDIT")
print("=" * 100)

print("Donors:", len(observed_donors))
print("Layers:", len(EXPECTED_LAYERS))
print(
    "Target cells excluded for missing layer annotation:",
    len(exclusions),
)
print(
    "Donor × layer combinations:",
    n_total,
)

print(
    "Passing min_ncells=10:",
    n_pass,
)

print(
    "Failing min_ncells=10:",
    n_fail,
)


if n_fail == 0:

    print(
        "\nFINAL STATUS: PASS — "
        "all 23 × 7 donor-layer combinations "
        "satisfy min_ncells=10"
    )

else:

    print(
        "\nFINAL STATUS: CONDITIONAL PASS — "
        f"{n_fail} donor-layer combinations "
        "have <10 cells and must be excluded"
    )


print("\nOutput directory:")
print(OUT)
