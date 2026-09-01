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


# ==============================================================================
# LOAD
# ==============================================================================

section("LOAD ANNOTATED XENIUM")

adata = sc.read_h5ad(
    H5AD,
    backed="r",
)

print("Shape:", adata.shape)

required = [
    "BrNum",
    "Dx",
    "cell_type_annotation",
]

missing = [
    x for x in required
    if x not in adata.obs.columns
]

if missing:
    raise RuntimeError(
        f"Missing obs columns: {missing}"
    )


obs = adata.obs[
    required
].copy()

obs["BrNum"] = obs["BrNum"].astype(str)
obs["Dx"] = obs["Dx"].astype(str)
obs["cell_type_annotation"] = (
    obs["cell_type_annotation"]
    .astype(str)
)


# ==============================================================================
# N23
# ==============================================================================

section("CANONICAL N23")

obs = obs[
    obs["BrNum"] != "Br6432"
].copy()

print(
    "Donors:",
    obs["BrNum"].nunique()
)

print(
    "Cells:",
    len(obs)
)

print(
    "Br6432 excluded:",
    "Br6432" not in set(obs["BrNum"])
)


donor_dx = (
    obs[
        ["BrNum", "Dx"]
    ]
    .drop_duplicates()
)

print("\nDiagnosis by donor:")
print(
    donor_dx["Dx"]
    .value_counts()
    .sort_index()
)


# ==============================================================================
# CELL TYPE LABELS
# ==============================================================================

section("CELL TYPE LABELS")

labels = sorted(
    obs["cell_type_annotation"]
    .dropna()
    .unique()
)

print(
    "Number of observed cell-type labels:",
    len(labels)
)

print("\nObserved labels:")

for x in labels:
    print(" -", x)


observed = set(labels)
expected = set(PAPER_CELL_TYPES)

missing_from_ours = sorted(
    expected - observed
)

extra_in_ours = sorted(
    observed - expected
)


print("\nExpected paper cell types:", len(expected))

print(
    "Missing paper labels:",
    missing_from_ours
)

print(
    "Extra labels in GeneBridge:",
    extra_in_ours
)

print(
    "Exact paper label set:",
    observed == expected
)


# ==============================================================================
# TOTAL CELL COUNTS
# ==============================================================================

section("TOTAL CELLS PER CELL TYPE")

total_counts = (
    obs["cell_type_annotation"]
    .value_counts()
    .rename_axis("cell_type")
    .reset_index(name="n_cells")
)

print(
    total_counts.to_string(
        index=False
    )
)

total_counts.to_csv(
    OUT / "09_celltype_total_cell_counts.csv",
    index=False,
)


# ==============================================================================
# DONOR × CELL TYPE COUNTS
# ==============================================================================

section("DONOR x CELL TYPE COUNTS")

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
    .reset_index(name="n_cells")
)

counts["passes_min10"] = (
    counts["n_cells"] >= 10
)

counts.to_csv(
    OUT / "09_donor_celltype_counts.csv",
    index=False,
)


# All donor × label combinations, including missing combinations.

all_donors = sorted(
    obs["BrNum"].unique()
)

all_celltypes = sorted(
    obs["cell_type_annotation"].unique()
)

grid = pd.MultiIndex.from_product(
    [
        all_donors,
        all_celltypes,
    ],
    names=[
        "BrNum",
        "cell_type",
    ],
).to_frame(index=False)

dx_map = (
    donor_dx
    .set_index("BrNum")["Dx"]
)

grid["Dx"] = grid["BrNum"].map(dx_map)

tmp = counts.rename(
    columns={
        "cell_type_annotation": "cell_type"
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

grid.to_csv(
    OUT / "09_all_donor_celltype_coverage.csv",
    index=False,
)


# ==============================================================================
# SUMMARY BY CELL TYPE
# ==============================================================================

section("PSEUDOBULK ELIGIBILITY BY CELL TYPE")

summary = (
    grid.groupby(
        "cell_type",
        observed=True,
    )
    .agg(
        donors_total=("BrNum", "nunique"),
        donors_with_cells=(
            "n_cells",
            lambda x: int((x > 0).sum())
        ),
        donors_passing_min10=(
            "passes_min10",
            "sum"
        ),
        min_cells=("n_cells", "min"),
        median_cells=("n_cells", "median"),
        max_cells=("n_cells", "max"),
    )
    .reset_index()
)


# SCZ/NTC donor coverage after min10

for dx in ["NTC", "SCZ"]:

    sub = grid[
        grid["Dx"] == dx
    ]

    dx_pass = (
        sub.groupby(
            "cell_type",
            observed=True,
        )["passes_min10"]
        .sum()
    )

    summary[
        f"{dx}_donors_passing_min10"
    ] = (
        summary["cell_type"]
        .map(dx_pass)
        .fillna(0)
        .astype(int)
    )


summary = summary.sort_values(
    "cell_type"
)

print(
    summary.to_string(
        index=False
    )
)

summary.to_csv(
    OUT / "09_celltype_pseudobulk_eligibility_summary.csv",
    index=False,
)


# ==============================================================================
# EXPECTED NUMBER OF PSEUDOBULKS
# ==============================================================================

section("EXPECTED N23 CELL-TYPE PSEUDOBULKS")

eligible = grid[
    grid["passes_min10"]
]

print(
    "Eligible donor-celltype pseudobulks:",
    len(eligible)
)

print(
    "Maximum possible with 23 donors x",
    len(all_celltypes),
    "cell types:",
    23 * len(all_celltypes),
)

print(
    "\nEligible PBs by diagnosis:"
)

print(
    eligible["Dx"]
    .value_counts()
    .sort_index()
)


# ==============================================================================
# PROBLEMATIC CELL TYPES
# ==============================================================================

section("CELL TYPES WITH INCOMPLETE DONOR COVERAGE")

problem = summary[
    summary["donors_passing_min10"] < 23
].copy()

if len(problem) == 0:
    print(
        "All cell types pass min_ncells=10 "
        "for all 23 donors."
    )
else:
    print(
        problem.to_string(
            index=False
        )
    )


# ==============================================================================
# FINAL
# ==============================================================================

section("FINAL STATUS")

print("Canonical donors: 23")
print("Paper expected cell types:", len(PAPER_CELL_TYPES))
print("Observed cell types:", len(labels))
print("Exact paper label set:", observed == expected)

print(
    "Eligible donor-celltype pseudobulks:",
    len(eligible)
)

print(
    "\nFINAL STATUS: N23 CELL-TYPE PSEUDOBULK AUDIT COMPLETE"
)
