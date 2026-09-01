#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MANIFEST = ROOT / "data/metadata/envi_production_23donors.tsv"
DONOR_META = ROOT / "data/metadata/xenium_DE_metadata_23.csv"

ANNOTATED = (
    ROOT
    / "data/processed/xenium"
    / "xenium_N24_layer_celltype_annotated.h5ad"
)

COVERAGE = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "celltype_N23"
    / "09c_N23_donor_celltype_coverage.csv"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "celltype_N23"
    / "pseudobulk"
    / "per_donor"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


CELL_TYPES = [
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

MIN_NCELLS = 10
CHUNK = 512


def section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def sum_block(x):
    if sparse.issparse(x):
        return (
            np.asarray(
                x.sum(axis=0)
            )
            .ravel()
            .astype(np.float64)
        )

    return np.asarray(
        x,
        dtype=np.float64,
    ).sum(
        axis=0,
        dtype=np.float64,
    )


# =============================================================================
# ARGUMENT
# =============================================================================

parser = argparse.ArgumentParser()

parser.add_argument(
    "--index",
    required=True,
    type=int,
    help="1-based row in the 23-donor ENVI manifest",
)

args = parser.parse_args()

if not 1 <= args.index <= 23:
    raise ValueError(
        "--index must be 1-23"
    )


# =============================================================================
# RESOLVE DONOR
# =============================================================================

section("RESOLVE DONOR")

manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
)

manifest["donor"] = (
    manifest["donor"]
    .astype(str)
    .str.strip()
)

manifest["Dx"] = (
    manifest["Dx"]
    .astype(str)
    .str.strip()
    .str.upper()
)

if len(manifest) != 23:
    raise RuntimeError(
        f"Expected 23 manifest donors; found {len(manifest)}"
    )


r = manifest.iloc[
    args.index - 1
]

donor = str(r["donor"])
dx = str(r["Dx"])
experiment = str(r["experiment"])
target_path = Path(r["target"])


print("Array index:", args.index)
print("Donor:", donor)
print("Dx:", dx)
print("Experiment:", experiment)
print("Target:", target_path)


# =============================================================================
# TRUSTED DONOR METADATA
# =============================================================================

section("CANONICAL DONOR METADATA")

dm = pd.read_csv(
    DONOR_META
)

dm["BrNum"] = (
    dm["BrNum"]
    .astype(str)
    .str.strip()
)

dm = dm.loc[
    dm["BrNum"] == donor
]

if len(dm) != 1:
    raise RuntimeError(
        f"{donor}: canonical metadata rows = {len(dm)}"
    )

dm = dm.iloc[0]


if str(dm["Dx"]).upper() != dx:
    raise RuntimeError(
        f"{donor}: diagnosis mismatch "
        f"manifest={dx}, metadata={dm['Dx']}"
    )


print("Dx:", dm["Dx"])
print("Age:", dm["Age"])
print("Sex:", dm["Sex"])
print("slide_id:", dm["slide_id"])
print("run_date:", dm["run_date"])


# =============================================================================
# EXPECTED CELL-TYPE COVERAGE
# =============================================================================

section("EXPECTED DONOR x CELL-TYPE COVERAGE")

coverage = pd.read_csv(
    COVERAGE
)

coverage["BrNum"] = (
    coverage["BrNum"]
    .astype(str)
)

expected = coverage.loc[
    coverage["BrNum"] == donor
].copy()


if len(expected) != 12:
    raise RuntimeError(
        f"{donor}: expected 12 coverage rows; "
        f"found {len(expected)}"
    )


if not expected["passes_min10"].all():
    raise RuntimeError(
        f"{donor}: one or more cell types failed min10 audit"
    )


expected_counts = dict(
    zip(
        expected["cell_type"],
        expected["n_cells"],
    )
)


if set(expected_counts) != set(CELL_TYPES):
    raise RuntimeError(
        f"{donor}: expected cell-type set mismatch"
    )


for ct in CELL_TYPES:
    print(
        f"{ct:20s} "
        f"{int(expected_counts[ct]):8,d}"
    )


# =============================================================================
# LOAD CELL-TYPE ANNOTATION
# =============================================================================

section("LOAD CELL-TYPE ANNOTATION")

ann = ad.read_h5ad(
    ANNOTATED,
    backed="r",
)

mask = (
    ann.obs["BrNum"]
    .astype(str)
    .eq(donor)
)

annotation = (
    ann.obs.loc[
        mask,
        [
            "BrNum",
            "cell_type_annotation",
        ],
    ]
    .copy()
)

ann.file.close()


annotation["BrNum"] = (
    annotation["BrNum"]
    .astype(str)
)

annotation["cell_type_annotation"] = (
    annotation["cell_type_annotation"]
    .astype(str)
)


print(
    "Annotated donor cells:",
    len(annotation)
)

print("\nAnnotation cell-type counts:")

print(
    annotation[
        "cell_type_annotation"
    ]
    .value_counts()
    .reindex(
        CELL_TYPES,
        fill_value=0,
    )
    .to_string()
)


unknown_labels = (
    set(
        annotation[
            "cell_type_annotation"
        ]
    )
    - set(CELL_TYPES)
)

if unknown_labels:
    raise RuntimeError(
        f"{donor}: unexpected cell-type labels: "
        f"{sorted(unknown_labels)}"
    )


# =============================================================================
# OPEN ORIGINAL TARGET + ENVI
# =============================================================================

section("OPEN TARGET + ENVI")

after_path = (
    ROOT
    / "data/processed/imputation_full"
    / experiment
    / "envi"
    / donor
    / f"spatial_data_xenium_{donor}_ENVI_full_transcriptome.h5ad"
)


if not target_path.exists():
    raise FileNotFoundError(
        target_path
    )

if not after_path.exists():
    raise FileNotFoundError(
        after_path
    )


target = ad.read_h5ad(
    target_path,
    backed="r",
)

after = ad.read_h5ad(
    after_path,
    backed="r",
)


print(
    "Target:",
    target.shape
)

print(
    "ENVI:",
    after.shape
)


target_ids = pd.Index(
    target.obs_names.astype(str)
)

after_ids = pd.Index(
    after.obs_names.astype(str)
)


if len(target_ids) != len(after_ids):
    raise RuntimeError(
        f"{donor}: target/ENVI cell counts differ"
    )

if not np.array_equal(
    target_ids.to_numpy(),
    after_ids.to_numpy(),
):
    raise RuntimeError(
        f"{donor}: target/ENVI cell IDs or order differ"
    )


ann_ids = pd.Index(
    annotation.index.astype(str)
)


target_only = target_ids.difference(
    ann_ids
)

annotated_only = ann_ids.difference(
    target_ids
)


print(
    "Target cells:",
    len(target_ids)
)

print(
    "Target-only unannotated:",
    len(target_only)
)

print(
    "Annotated-only:",
    len(annotated_only)
)


if len(annotated_only) > 0:
    raise RuntimeError(
        f"{donor}: {len(annotated_only)} annotated cells "
        "are absent from ENVI target"
    )


# =============================================================================
# GENE VALIDATION
# =============================================================================

section("GENE VALIDATION")

if target.n_vars != 300:
    raise RuntimeError(
        f"{donor}: target expected 300 genes; "
        f"found {target.n_vars}"
    )

if after.n_vars != 34987:
    raise RuntimeError(
        f"{donor}: ENVI expected 34,987 genes; "
        f"found {after.n_vars}"
    )


original_genes = (
    target.var_names
    .astype(str)
    .tolist()
)

full_genes = (
    after.var_names
    .astype(str)
    .tolist()
)


if len(set(original_genes)) != 300:
    raise RuntimeError(
        f"{donor}: duplicate target genes"
    )

if len(set(full_genes)) != 34987:
    raise RuntimeError(
        f"{donor}: duplicate ENVI genes"
    )


if "expression_source" not in after.var.columns:
    raise RuntimeError(
        f"{donor}: expression_source missing"
    )


source = (
    after.var[
        "expression_source"
    ]
    .astype(str)
)

print("\nENVI expression source:")

print(
    source.value_counts()
)


if (
    (source == "measured_xenium").sum()
    != 300
):
    raise RuntimeError(
        f"{donor}: expected 300 measured genes"
    )


if (
    (source == "envi_imputed").sum()
    != 34687
):
    raise RuntimeError(
        f"{donor}: expected 34,687 imputed genes"
    )


full_gene_index = pd.Index(
    full_genes
)

measured_idx = (
    full_gene_index
    .get_indexer(
        original_genes
    )
)


if np.any(
    measured_idx < 0
):
    missing = np.asarray(
        original_genes
    )[
        measured_idx < 0
    ]

    raise RuntimeError(
        f"{donor}: measured genes missing from ENVI: "
        f"{missing[:20]}"
    )


if not all(
    source.iloc[j]
    == "measured_xenium"
    for j in measured_idx
):
    raise RuntimeError(
        f"{donor}: some target genes are not "
        "marked measured_xenium"
    )


print(
    "Original Xenium genes:",
    len(original_genes)
)

print(
    "ENVI full genes:",
    len(full_genes)
)

print(
    "Original genes present as measured:",
    len(measured_idx)
)


# =============================================================================
# MATRIX SOURCES
# =============================================================================

section("MATRIX SOURCES")


# Original Xenium measured raw counts.
if "counts" in target.layers:
    original_matrix = target.layers[
        "counts"
    ]
    print(
        "Original Xenium source: target.layers['counts']"
    )
else:
    original_matrix = target.X
    print(
        "WARNING: target.layers['counts'] absent; "
        "using target.X"
    )


# Preserved measured Xenium values in ENVI production object.
measured_matrix = after.X

print(
    "ENVI measured source: after.X"
)


# Full ENVI transcriptome on count scale.
if "count_scale" not in after.layers:
    raise RuntimeError(
        f"{donor}: ENVI count_scale layer missing"
    )

full_matrix = after.layers[
    "count_scale"
]

print(
    "ENVI full source: after.layers['count_scale']"
)


# =============================================================================
# ALLOCATE CELL-TYPE TOTALS
# =============================================================================

section("ALLOCATE CELL-TYPE TOTALS")


original_totals = {
    ct: np.zeros(
        300,
        dtype=np.float64,
    )
    for ct in CELL_TYPES
}


measured_totals = {
    ct: np.zeros(
        300,
        dtype=np.float64,
    )
    for ct in CELL_TYPES
}


full_totals = {
    ct: np.zeros(
        34987,
        dtype=np.float64,
    )
    for ct in CELL_TYPES
}


celltype_cells = {
    ct: 0
    for ct in CELL_TYPES
}


matched_original_total = np.zeros(
    300,
    dtype=np.float64,
)

matched_measured_total = np.zeros(
    300,
    dtype=np.float64,
)

matched_full_total = np.zeros(
    34987,
    dtype=np.float64,
)


# =============================================================================
# CHUNK THROUGH CELLS
# =============================================================================

section("AGGREGATE CELLS")

for start in range(
    0,
    after.n_obs,
    CHUNK,
):

    stop = min(
        start + CHUNK,
        after.n_obs,
    )


    chunk_ids = target_ids[
        start:stop
    ]


    keep = chunk_ids.isin(
        ann_ids
    )


    if not keep.any():
        continue


    keep_positions = np.flatnonzero(
        keep
    )


    kept_ids = chunk_ids[
        keep
    ]


    kept_ann = annotation.loc[
        kept_ids,
        [
            "cell_type_annotation",
        ],
    ]


    # -------------------------------------------------------------------------
    # ORIGINAL 300
    # -------------------------------------------------------------------------

    o_all = original_matrix[
        start:stop,
        :
    ]

    if sparse.issparse(
        o_all
    ):
        o = o_all[
            keep_positions,
            :
        ]
    else:
        o = np.asarray(
            o_all,
            dtype=np.float64,
        )[
            keep_positions,
            :
        ]


    # -------------------------------------------------------------------------
    # ENVI MEASURED 300
    # -------------------------------------------------------------------------

    m_all = measured_matrix[
        start:stop,
        measured_idx,
    ]


    if sparse.issparse(
        m_all
    ):
        m = m_all[
            keep_positions,
            :
        ]
    else:
        m = np.asarray(
            m_all,
            dtype=np.float64,
        )[
            keep_positions,
            :
        ]


    # -------------------------------------------------------------------------
    # ENVI FULL 34,987
    # -------------------------------------------------------------------------

    f_all = full_matrix[
        start:stop,
        :
    ]


    if sparse.issparse(
        f_all
    ):
        f = f_all[
            keep_positions,
            :
        ]
    else:
        f = np.asarray(
            f_all,
            dtype=np.float64,
        )[
            keep_positions,
            :
        ]


    matched_original_total += sum_block(
        o
    )

    matched_measured_total += sum_block(
        m
    )

    matched_full_total += sum_block(
        f
    )


    # -------------------------------------------------------------------------
    # SPLIT BY CELL TYPE
    # -------------------------------------------------------------------------

    for ct in CELL_TYPES:

        cm = (
            kept_ann[
                "cell_type_annotation"
            ]
            .eq(ct)
            .to_numpy()
        )

        if not cm.any():
            continue


        n = int(
            cm.sum()
        )


        celltype_cells[
            ct
        ] += n


        original_totals[
            ct
        ] += sum_block(
            o[
                cm,
                :
            ]
        )


        measured_totals[
            ct
        ] += sum_block(
            m[
                cm,
                :
            ]
        )


        full_totals[
            ct
        ] += sum_block(
            f[
                cm,
                :
            ]
        )


    if (
        start == 0
        or stop == after.n_obs
        or start % (CHUNK * 20) == 0
    ):

        print(
            f"cells {stop:,}/{after.n_obs:,}",
            flush=True,
        )


# =============================================================================
# CELL COUNT VALIDATION
# =============================================================================

section("CELL COUNT VALIDATION")


for ct in CELL_TYPES:

    observed = int(
        celltype_cells[ct]
    )

    expected_n = int(
        expected_counts[ct]
    )


    print(
        f"{ct:20s} "
        f"observed={observed:8,d} "
        f"expected={expected_n:8,d}"
    )


    if observed != expected_n:
        raise RuntimeError(
            f"{donor} {ct}: "
            f"observed cells={observed}, "
            f"expected={expected_n}"
        )


    if observed < MIN_NCELLS:
        raise RuntimeError(
            f"{donor} {ct}: "
            f"{observed} cells < {MIN_NCELLS}"
        )


print(
    "\nAll 12 cell-type counts match audit: TRUE"
)


# =============================================================================
# SUM-CONSERVATION VALIDATION
# =============================================================================

section("SUM CONSERVATION")


original_pb_sum = (
    np.stack(
        [
            original_totals[ct]
            for ct in CELL_TYPES
        ]
    )
    .sum(axis=0)
)


measured_pb_sum = (
    np.stack(
        [
            measured_totals[ct]
            for ct in CELL_TYPES
        ]
    )
    .sum(axis=0)
)


full_pb_sum = (
    np.stack(
        [
            full_totals[ct]
            for ct in CELL_TYPES
        ]
    )
    .sum(axis=0)
)


if not np.allclose(
    original_pb_sum,
    matched_original_total,
    rtol=1e-10,
    atol=1e-6,
):
    raise RuntimeError(
        f"{donor}: original PB sum mismatch"
    )


if not np.allclose(
    measured_pb_sum,
    matched_measured_total,
    rtol=1e-10,
    atol=1e-6,
):
    raise RuntimeError(
        f"{donor}: measured PB sum mismatch"
    )


if not np.allclose(
    full_pb_sum,
    matched_full_total,
    rtol=1e-10,
    atol=1e-6,
):
    raise RuntimeError(
        f"{donor}: full PB sum mismatch"
    )


print(
    "Original sum conservation: TRUE"
)

print(
    "ENVI measured sum conservation: TRUE"
)

print(
    "ENVI full sum conservation: TRUE"
)


# =============================================================================
# BUILD OUTPUT TABLES
# =============================================================================

section("BUILD OUTPUT TABLES")


original_rows = []
measured_rows = []
full_rows = []
metadata_rows = []


for ct in CELL_TYPES:

    safe_ct = (
        ct.replace(
            "/",
            "-"
        )
        .replace(
            ":",
            ""
        )
        .replace(
            ",",
            ""
        )
        .replace(
            " ",
            "_"
        )
    )


    pb_id = (
        f"{donor}::{ct}"
    )


    ov = original_totals[
        ct
    ]

    mv = measured_totals[
        ct
    ]

    fv = full_totals[
        ct
    ]


    for label, vec in [
        ("original", ov),
        ("measured", mv),
        ("full", fv),
    ]:

        if not np.isfinite(
            vec
        ).all():

            raise RuntimeError(
                f"{pb_id}: "
                f"{label} contains non-finite values"
            )


        if np.any(
            vec < 0
        ):

            raise RuntimeError(
                f"{pb_id}: "
                f"{label} contains negative values"
            )


    original_rows.append(
        pd.Series(
            ov,
            index=original_genes,
            name=pb_id,
        )
    )


    measured_rows.append(
        pd.Series(
            mv,
            index=original_genes,
            name=pb_id,
        )
    )


    full_rows.append(
        pd.Series(
            fv,
            index=full_genes,
            name=pb_id,
        )
    )


    metadata_rows.append(
        {
            "pseudobulk_id": pb_id,
            "BrNum": donor,
            "Dx": dx,
            "Age": float(
                dm["Age"]
            ),
            "Sex": str(
                dm["Sex"]
            ),
            "slide_id": str(
                dm["slide_id"]
            ),
            "run_date": str(
                dm["run_date"]
            ),
            "cell_type": ct,
            "annots": ct,
            "n_cells": int(
                celltype_cells[ct]
            ),
            "original300_library":
                float(
                    ov.sum()
                ),
            "envi_measured300_library":
                float(
                    mv.sum()
                ),
            "envi_full34987_library":
                float(
                    fv.sum()
                ),
        }
    )


original_df = pd.DataFrame(
    original_rows
)

measured_df = pd.DataFrame(
    measured_rows
)

full_df = pd.DataFrame(
    full_rows
)

metadata_df = pd.DataFrame(
    metadata_rows
)


# =============================================================================
# ORIGINAL VS ENVI-MEASURED VALIDATION
# =============================================================================

section("ORIGINAL VS ENVI MEASURED CONTROL")


print(
    "Shapes:",
    original_df.shape,
    measured_df.shape,
)


same_exact = np.array_equal(
    original_df.to_numpy(),
    measured_df.to_numpy(),
)


same_close = np.allclose(
    original_df.to_numpy(),
    measured_df.to_numpy(),
    rtol=0,
    atol=1e-8,
)


max_abs_diff = float(
    np.max(
        np.abs(
            original_df.to_numpy()
            -
            measured_df.to_numpy()
        )
    )
)


print(
    "Original vs ENVI measured exact:",
    same_exact
)

print(
    "Original vs ENVI measured allclose:",
    same_close
)

print(
    "Maximum absolute difference:",
    max_abs_diff
)


# =============================================================================
# SAVE
# =============================================================================

section("SAVE DONOR OUTPUTS")


prefix = (
    OUT
    / donor
)


original_file = Path(
    str(prefix)
    + "_original300_celltype_pseudobulk.csv.gz"
)

measured_file = Path(
    str(prefix)
    + "_ENVI_measured300_celltype_pseudobulk.csv.gz"
)

full_file = Path(
    str(prefix)
    + "_ENVI_full34987_celltype_pseudobulk.csv.gz"
)

meta_file = Path(
    str(prefix)
    + "_celltype_pseudobulk_metadata.csv"
)

exclusion_file = Path(
    str(prefix)
    + "_excluded_target_cells.csv"
)


original_df.to_csv(
    original_file,
    compression="gzip",
)

measured_df.to_csv(
    measured_file,
    compression="gzip",
)

full_df.to_csv(
    full_file,
    compression="gzip",
)

metadata_df.to_csv(
    meta_file,
    index=False,
)


pd.DataFrame(
    {
        "BrNum": donor,
        "cell_id": [
            str(x)
            for x in target_only
        ],
        "reason":
            "missing_celltype_annotation",
    }
).to_csv(
    exclusion_file,
    index=False,
)


# =============================================================================
# CLOSE FILES
# =============================================================================

target.file.close()
after.file.close()


# =============================================================================
# FINAL
# =============================================================================

section("FINAL STATUS")


print("Donor:", donor)
print("Diagnosis:", dx)
print("Cell types:", len(CELL_TYPES))
print("Pseudobulks:", len(metadata_df))

print(
    "Original matrix:",
    original_df.shape
)

print(
    "ENVI measured matrix:",
    measured_df.shape
)

print(
    "ENVI full matrix:",
    full_df.shape
)

print(
    "Target-only excluded cells:",
    len(target_only)
)

print(
    "Original vs ENVI measured exact:",
    same_exact
)

print(
    "Original vs ENVI measured allclose:",
    same_close
)

print("\nSaved:")
print(original_file)
print(measured_file)
print(full_file)
print(meta_file)
print(exclusion_file)

print(
    "\nFINAL STATUS: "
    "DONOR CELL-TYPE PSEUDOBULK COMPLETE"
)
