#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

PER = (
    ROOT
    / "outputs/AIVC_beta/presplit_audit/per_donor"
)

OUT = (
    ROOT
    / "outputs/AIVC_beta/presplit_audit"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


CELL_TYPES = [
    "Oligo",
    "Endo",
    "Ast",
    "L2/3 Ex",
    "Mic",
    "In: SST, PVALB",
    "L6 Ex",
    "Ambig/In/Endo",
    "L4/5 Ex",
    "In: VIP, LAMP5",
    "Ambig/Oligo",
    "L5 Ex",
]

LAYERS = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WMtz",
    "WM",
]


def read_all(pattern):

    files = sorted(
        PER.glob(pattern)
    )

    if not files:
        return pd.DataFrame()

    return pd.concat(
        [
            pd.read_csv(f)
            for f in files
        ],
        ignore_index=True,
    )


donors = read_all(
    "Br*_donor_audit.csv"
)

celltypes = read_all(
    "Br*_celltype_counts.csv"
)

layers = read_all(
    "Br*_layer_counts.csv"
)

expr = read_all(
    "Br*_expression_qc.csv"
)


print("=" * 100)
print("AIVC PRE-SPLIT AUDIT AGGREGATION")
print("=" * 100)

print(
    "Per-donor audit files:",
    len(donors),
)


if len(donors) != 24:

    raise RuntimeError(
        f"Expected 24 completed donors; found {len(donors)}"
    )


donors = donors.sort_values(
    [
        "Dx",
        "BrNum",
    ]
)

celltypes = celltypes.sort_values(
    [
        "Dx",
        "BrNum",
    ]
)

layers = layers.sort_values(
    [
        "Dx",
        "BrNum",
    ]
)


# ==============================================================================
# CELL-TYPE COVERAGE
# ==============================================================================

coverage_rows = []

for ct in CELL_TYPES:

    ntc = celltypes[
        celltypes["Dx"] == "NTC"
    ]

    scz = celltypes[
        celltypes["Dx"] == "SCZ"
    ]

    positive = (
        celltypes[ct] > 0
    )

    coverage_rows.append(
        {
            "cell_type": ct,

            "total_cells":
                int(
                    celltypes[
                        ct
                    ].sum()
                ),

            "donors_present":
                int(
                    positive.sum()
                ),

            "NTC_donors_present":
                int(
                    (
                        ntc[
                            ct
                        ] > 0
                    ).sum()
                ),

            "SCZ_donors_present":
                int(
                    (
                        scz[
                            ct
                        ] > 0
                    ).sum()
                ),

            "NTC_cells":
                int(
                    ntc[
                        ct
                    ].sum()
                ),

            "SCZ_cells":
                int(
                    scz[
                        ct
                    ].sum()
                ),

            "min_cells_when_present":
                int(
                    celltypes.loc[
                        positive,
                        ct,
                    ].min()
                )
                if positive.any()
                else 0,
        }
    )

coverage = pd.DataFrame(
    coverage_rows
)


# ==============================================================================
# LAYER COVERAGE
# ==============================================================================

layer_rows = []

for layer in LAYERS:

    ntc = layers[
        layers["Dx"] == "NTC"
    ]

    scz = layers[
        layers["Dx"] == "SCZ"
    ]

    layer_rows.append(
        {
            "layer": layer,

            "total_cells":
                int(
                    layers[
                        layer
                    ].sum()
                ),

            "donors_present":
                int(
                    (
                        layers[
                            layer
                        ] > 0
                    ).sum()
                ),

            "NTC_donors_present":
                int(
                    (
                        ntc[
                            layer
                        ] > 0
                    ).sum()
                ),

            "SCZ_donors_present":
                int(
                    (
                        scz[
                            layer
                        ] > 0
                    ).sum()
                ),
        }
    )

layer_coverage = pd.DataFrame(
    layer_rows
)


# ==============================================================================
# TECHNICAL × DIAGNOSIS
# ==============================================================================

technical_rows = []

for col in [
    "slide_id",
    "run_date",
    "CaptureArea",
]:

    if col not in donors.columns:
        continue

    ntc_levels = set(
        donors.loc[
            donors["Dx"] == "NTC",
            col,
        ]
        .dropna()
        .astype(str)
    )

    scz_levels = set(
        donors.loc[
            donors["Dx"] == "SCZ",
            col,
        ]
        .dropna()
        .astype(str)
    )

    shared = (
        ntc_levels
        & scz_levels
    )

    technical_rows.append(
        {
            "variable": col,
            "NTC_levels":
                len(
                    ntc_levels
                ),
            "SCZ_levels":
                len(
                    scz_levels
                ),
            "shared_levels":
                len(
                    shared
                ),
            "shared_values":
                ";".join(
                    sorted(
                        shared
                    )
                ),
        }
    )

technical = pd.DataFrame(
    technical_rows
)


# ==============================================================================
# STATUS CHECKS
# ==============================================================================

checks = []


def add(
    name,
    passed,
    detail,
    warning=False,
):

    if passed:
        status = "PASS"
    elif warning:
        status = "WARN"
    else:
        status = "FAIL"

    checks.append(
        {
            "check": name,
            "status": status,
            "detail": detail,
        }
    )


dx_counts = (
    donors[
        "Dx"
    ]
    .value_counts()
    .to_dict()
)

add(
    "24 ENVI donors",
    len(donors) == 24,
    f"found={len(donors)}",
)

add(
    "12 NTC + 12 SCZ",
    (
        dx_counts.get(
            "NTC",
            0,
        ) == 12
        and
        dx_counts.get(
            "SCZ",
            0,
        ) == 12
    ),
    str(dx_counts),
)

add(
    "Donor identity per H5AD",
    donors[
        "donor_matches_path"
    ].astype(bool).all(),
    (
        f"{int(donors['donor_matches_path'].astype(bool).sum())}"
        f"/{len(donors)}"
    ),
)

add(
    "Exact ENVI-annotation join",
    np.allclose(
        donors[
            "join_pct"
        ],
        100.0,
    ),
    (
        f"minimum="
        f"{donors['join_pct'].min():.4f}%"
    ),
)

add(
    "Unique ENVI cell IDs",
    (
        donors[
            "duplicate_envi_ids"
        ].sum()
        == 0
    ),
    (
        f"duplicates="
        f"{int(donors['duplicate_envi_ids'].sum())}"
    ),
)

add(
    "Spatial coordinates finite",
    donors[
        "spatial_finite"
    ].astype(bool).all(),
    (
        f"{int(donors['spatial_finite'].astype(bool).sum())}"
        f"/24"
    ),
)

gene_hashes = (
    donors[
        "gene_hash"
    ]
    .astype(str)
    .nunique()
)

add(
    "Same gene order",
    gene_hashes == 1,
    f"unique_gene_hashes={gene_hashes}",
)

ref_overlap = (
    donors.loc[
        donors[
            "reference_donor_overlap"
        ].astype(bool),
        "BrNum",
    ]
    .tolist()
)

add(
    "snRNA donor overlap",
    len(
        ref_overlap
    ) == 0,
    (
        "none"
        if len(
            ref_overlap
        ) == 0
        else ",".join(
            ref_overlap
        )
    ),
    warning=True,
)

if not expr.empty:

    add(
        "Expression sampled values finite",
        expr[
            "all_finite"
        ]
        .fillna(
            False
        )
        .astype(bool)
        .all(),
        (
            f"{int(expr['all_finite'].fillna(False).astype(bool).sum())}"
            f"/{len(expr)}"
        ),
        warning=True,
    )


status = pd.DataFrame(
    checks
)


# ==============================================================================
# SAVE
# ==============================================================================

donors.to_csv(
    OUT
    / "00_donor_audit.csv",
    index=False,
)

celltypes.to_csv(
    OUT
    / "00_celltype_counts_by_donor.csv",
    index=False,
)

coverage.to_csv(
    OUT
    / "00_celltype_donor_coverage.csv",
    index=False,
)

layers.to_csv(
    OUT
    / "00_layer_counts_by_donor.csv",
    index=False,
)

layer_coverage.to_csv(
    OUT
    / "00_layer_donor_coverage.csv",
    index=False,
)

expr.to_csv(
    OUT
    / "00_expression_sample_qc.csv",
    index=False,
)

technical.to_csv(
    OUT
    / "00_technical_dx_overlap.csv",
    index=False,
)

status.to_csv(
    OUT
    / "00_presplit_status.csv",
    index=False,
)


summary = []

summary.append(
    "AIVC BETA PRE-SPLIT AUDIT"
)

summary.append(
    "=" * 100
)

summary.append("")

summary.append(
    f"Donors: {len(donors)}"
)

summary.append(
    f"NTC: {(donors['Dx'] == 'NTC').sum()}"
)

summary.append(
    f"SCZ: {(donors['Dx'] == 'SCZ').sum()}"
)

summary.append(
    f"Total cells: {int(donors['n_cells'].sum())}"
)

summary.append("")

summary.append(
    "CHECKS"
)

summary.append(
    "-" * 100
)

for _, row in status.iterrows():

    summary.append(
        f"{row['status']:5s} | "
        f"{row['check']:<35s} | "
        f"{row['detail']}"
    )

summary.append("")

summary.append(
    "REFERENCE OVERLAP"
)

summary.append(
    "-" * 100
)

summary.append(
    ", ".join(
        ref_overlap
    )
    if ref_overlap
    else "None"
)

summary.append("")

summary.append(
    "CELL-TYPE COVERAGE"
)

summary.append(
    "-" * 100
)

summary.append(
    coverage.to_string(
        index=False
    )
)

summary.append("")

summary.append(
    "TECHNICAL × DX"
)

summary.append(
    "-" * 100
)

summary.append(
    technical.to_string(
        index=False
    )
)

summary_file = (
    OUT
    / "00_presplit_summary.txt"
)

summary_file.write_text(
    "\n".join(
        summary
    )
    + "\n"
)

print(
    summary_file.read_text()
)

