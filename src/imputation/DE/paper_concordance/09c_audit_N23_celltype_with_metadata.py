#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import scanpy as sc


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

H5AD = (
    ROOT
    / "data/processed/xenium"
    / "xenium_N24_layer_celltype_annotated.h5ad"
)

META = (
    ROOT
    / "data/metadata"
    / "xenium_DE_metadata_23.csv"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "celltype_N23"
)

OUT.mkdir(parents=True, exist_ok=True)


PAPER_CELL_TYPES = [
    "Oligo",
    "Ambig/Oligo",
    "Mic",
    "Ambig/In/Endo",
    "L5 Ex",
    "L6 Ex",
    "L4/5 Ex",
    "L2/3 Ex",
    "Ast",
    "Endo",
    "In: VIP, LAMP5",
    "In: SST, PVALB",
]


def section(x):
    print("\n" + "=" * 100)
    print(x)
    print("=" * 100)


# =============================================================================
# LOAD
# =============================================================================

section("LOAD INPUTS")

adata = sc.read_h5ad(
    H5AD,
    backed="r",
)

meta = pd.read_csv(META)

print("H5AD:", adata.shape)
print("Metadata:", meta.shape)


required_meta = [
    "BrNum",
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "run_date",
]

missing = [
    c for c in required_meta
    if c not in meta.columns
]

if missing:
    raise RuntimeError(
        f"Missing metadata columns: {missing}"
    )


# =============================================================================
# CANONICAL N23
# =============================================================================

section("CANONICAL N23")

obs = adata.obs[
    [
        "BrNum",
        "cell_type_annotation",
    ]
].copy()

obs["BrNum"] = obs["BrNum"].astype(str)

obs["cell_type_annotation"] = (
    obs["cell_type_annotation"]
    .astype(str)
)

meta["BrNum"] = meta["BrNum"].astype(str)


canonical_donors = set(
    meta["BrNum"]
)

obs = obs[
    obs["BrNum"].isin(
        canonical_donors
    )
].copy()


print(
    "Cells:",
    len(obs)
)

print(
    "Donors:",
    obs["BrNum"].nunique()
)

print(
    "Br6432 excluded:",
    "Br6432" not in set(obs["BrNum"])
)


if obs["BrNum"].nunique() != 23:
    raise RuntimeError(
        "Expected exactly 23 donors."
    )


if set(obs["BrNum"]) != canonical_donors:
    raise RuntimeError(
        "H5AD N23 donors do not exactly match metadata donors."
    )


# =============================================================================
# JOIN CANONICAL METADATA
# =============================================================================

section("JOIN CANONICAL METADATA")

obs = (
    obs.reset_index()
    .merge(
        meta[
            required_meta
        ],
        on="BrNum",
        how="left",
        validate="many_to_one",
    )
)

if obs["Dx"].isna().any():
    raise RuntimeError(
        "Missing Dx after metadata join."
    )

if obs["Age"].isna().any():
    raise RuntimeError(
        "Missing Age after metadata join."
    )

if obs["Sex"].isna().any():
    raise RuntimeError(
        "Missing Sex after metadata join."
    )


donor_meta = (
    obs[
        required_meta
    ]
    .drop_duplicates()
    .sort_values("BrNum")
)


print(
    "Unique donor metadata rows:",
    len(donor_meta)
)

print("\nDiagnosis by donor:")

print(
    donor_meta["Dx"]
    .value_counts()
    .sort_index()
)


if (
    donor_meta["Dx"]
    .value_counts()
    .to_dict()
    != {"NTC": 11, "SCZ": 12}
):
    raise RuntimeError(
        "Expected 11 NTC and 12 SCZ donors."
    )


# =============================================================================
# EXACT PAPER LABEL AUDIT
# =============================================================================

section("CELL-TYPE LABEL AUDIT")

observed_labels = set(
    obs["cell_type_annotation"]
    .unique()
)

expected_labels = set(
    PAPER_CELL_TYPES
)

print(
    "Observed cell types:",
    len(observed_labels)
)

print(
    "Expected paper cell types:",
    len(expected_labels)
)

print(
    "Missing paper labels:",
    sorted(expected_labels - observed_labels)
)

print(
    "Extra labels:",
    sorted(observed_labels - expected_labels)
)

print(
    "Exact paper label set:",
    observed_labels == expected_labels
)


if observed_labels != expected_labels:
    raise RuntimeError(
        "Cell-type labels do not exactly match paper."
    )


print("\nCell counts:")

print(
    obs["cell_type_annotation"]
    .value_counts()
    .sort_index()
    .to_string()
)


# =============================================================================
# DONOR x CELL TYPE
# =============================================================================

section("DONOR x CELL-TYPE COVERAGE")

counts = (
    obs.groupby(
        [
            "BrNum",
            "Dx",
            "cell_type_annotation",
        ],
        observed=True,
    )
    .size()
    .reset_index(
        name="n_cells"
    )
)


# complete donor × cell-type grid
grid = pd.MultiIndex.from_product(
    [
        sorted(canonical_donors),
        PAPER_CELL_TYPES,
    ],
    names=[
        "BrNum",
        "cell_type",
    ],
).to_frame(index=False)


donor_lookup = (
    meta.set_index("BrNum")
)


for c in [
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "run_date",
]:

    grid[c] = (
        grid["BrNum"]
        .map(donor_lookup[c])
    )


tmp = counts.rename(
    columns={
        "cell_type_annotation":
        "cell_type"
    }
)


grid = grid.merge(
    tmp[
        [
            "BrNum",
            "cell_type",
            "n_cells",
        ]
    ],
    on=[
        "BrNum",
        "cell_type",
    ],
    how="left",
)


grid["n_cells"] = (
    grid["n_cells"]
    .fillna(0)
    .astype(int)
)


grid["passes_min10"] = (
    grid["n_cells"] >= 10
)


print(
    "Possible donor-celltype PBs:",
    len(grid)
)

print(
    "Eligible donor-celltype PBs:",
    int(grid["passes_min10"].sum())
)


print("\nEligible PBs by diagnosis:")

print(
    grid.loc[
        grid["passes_min10"],
        "Dx"
    ]
    .value_counts()
    .sort_index()
)


# Expected:
# 11 NTC × 12 cell types = 132
# 12 SCZ × 12 cell types = 144

eligible_dx = (
    grid.loc[
        grid["passes_min10"],
        "Dx"
    ]
    .value_counts()
    .to_dict()
)


if eligible_dx != {
    "NTC": 132,
    "SCZ": 144,
}:

    raise RuntimeError(
        f"Unexpected PB diagnosis counts: {eligible_dx}"
    )


if grid["passes_min10"].sum() != 276:

    raise RuntimeError(
        "Expected all 276 donor-celltype PBs to pass min10."
    )


# =============================================================================
# CELL-TYPE SUMMARY
# =============================================================================

section("CELL-TYPE PSEUDOBULK SUMMARY")

summary = (
    grid.groupby(
        "cell_type",
        observed=True,
    )
    .agg(
        donors_total=(
            "BrNum",
            "nunique",
        ),
        donors_passing_min10=(
            "passes_min10",
            "sum",
        ),
        min_cells=(
            "n_cells",
            "min",
        ),
        median_cells=(
            "n_cells",
            "median",
        ),
        max_cells=(
            "n_cells",
            "max",
        ),
    )
    .reset_index()
)


for dx in [
    "NTC",
    "SCZ",
]:

    temp = (
        grid[
            grid["Dx"] == dx
        ]
        .groupby(
            "cell_type",
            observed=True,
        )["passes_min10"]
        .sum()
    )

    summary[
        f"{dx}_donors_passing_min10"
    ] = (
        summary["cell_type"]
        .map(temp)
        .astype(int)
    )


print(
    summary.to_string(
        index=False
    )
)


# =============================================================================
# SAVE
# =============================================================================

section("SAVE OUTPUTS")

grid.to_csv(
    OUT
    / "09c_N23_donor_celltype_coverage.csv",
    index=False,
)

summary.to_csv(
    OUT
    / "09c_N23_celltype_pseudobulk_summary.csv",
    index=False,
)

donor_meta.to_csv(
    OUT
    / "09c_N23_canonical_donor_metadata.csv",
    index=False,
)


# =============================================================================
# FINAL
# =============================================================================

section("FINAL STATUS")

print("Donors: 23")
print("NTC donors: 11")
print("SCZ donors: 12")
print("Cell types: 12")
print("Exact paper label set: TRUE")
print("Possible PBs: 276")
print("Eligible PBs: 276")
print("Minimum cells per PB: >= 10")

print(
    "\nFINAL STATUS: "
    "N23 CELL-TYPE AUDIT PASS"
)
