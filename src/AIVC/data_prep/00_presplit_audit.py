#!/usr/bin/env python3

from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd


# ======================================================================================
# PATHS
# ======================================================================================

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

ENVI_DIRS = [
    ROOT / "data/processed/imputation_full/ex1_ntc/envi",
    ROOT / "data/processed/imputation_full/ex2_scz/envi",
]

ANNOTATED = (
    ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad"
)

SNRNA_REFERENCE = (
    ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated/"
      "huuki_snrna_reference_full_allgenes.h5ad"
)

OUT = ROOT / "outputs/AIVC_beta/presplit_audit"
OUT.mkdir(parents=True, exist_ok=True)


EXPECTED_CELL_TYPES = [
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

EXPECTED_LAYERS = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WMtz",
    "WM",
]


# ======================================================================================
# HELPERS
# ======================================================================================

def section(title):
    print("\n" + "=" * 120)
    print(title)
    print("=" * 120)


def normalize_string_series(x):
    return x.astype(str).str.strip()


def unique_nonempty(x):
    vals = (
        normalize_string_series(x)
        .replace({"nan": np.nan, "None": np.nan, "": np.nan})
        .dropna()
        .unique()
        .tolist()
    )
    return sorted(vals)


def collapse_metadata(x):
    vals = unique_nonempty(x)

    if not vals:
        return np.nan

    return "|".join(vals)


def donor_from_path(path):
    # First preference: directory directly above file
    if re.fullmatch(r"Br\d+", path.parent.name):
        return path.parent.name

    m = re.search(r"(Br\d+)", path.name)

    if m:
        return m.group(1)

    raise RuntimeError(f"Could not infer donor from {path}")


# ======================================================================================
# FIND ENVI FILES
# ======================================================================================

section("LOCATE ENVI FULL-TRANSCRIPTOME FILES")

files = []

for d in ENVI_DIRS:

    if not d.exists():
        print("WARNING: directory does not exist:", d)
        continue

    files.extend(
        sorted(
            d.glob(
                "Br*/spatial_data_xenium_*_ENVI_full_transcriptome.h5ad"
            )
        )
    )

files = sorted(set(files))

print("Files found:", len(files))

for f in files:
    print(f)

if len(files) != 24:
    print(
        f"\nWARNING: expected 24 ENVI donor files, "
        f"but found {len(files)}."
    )


# ======================================================================================
# LOAD ANNOTATED XENIUM METADATA ONLY
# ======================================================================================

section("LOAD ANNOTATED XENIUM METADATA")

ann = ad.read_h5ad(
    ANNOTATED,
    backed="r",
)

ann_obs = ann.obs.copy()

if "BrNum" not in ann_obs.columns:
    raise RuntimeError("BrNum missing from annotated Xenium.")

ann_obs["BrNum"] = normalize_string_series(
    ann_obs["BrNum"]
)

print("Annotated Xenium cells:", len(ann_obs))
print(
    "Annotated Xenium donors:",
    ann_obs["BrNum"].nunique(),
)

print("\nAvailable annotation columns:")
print(ann_obs.columns.tolist())


# ======================================================================================
# REFERENCE DONORS
# ======================================================================================

section("LOAD snRNA REFERENCE DONOR IDS")

ref = ad.read_h5ad(
    SNRNA_REFERENCE,
    backed="r",
)

if "BrNum" not in ref.obs.columns:
    raise RuntimeError(
        "BrNum missing from snRNA reference."
    )

reference_donors = set(
    normalize_string_series(
        ref.obs["BrNum"]
    )
)

print("Reference donors:", len(reference_donors))
print(sorted(reference_donors))


# ======================================================================================
# AUDIT EACH ENVI DONOR
# ======================================================================================

section("PER-DONOR ENVI AUDIT")

donor_rows = []
gene_rows = []
expression_rows = []

celltype_count_records = []
layer_count_records = []

canonical_gene_names = None
canonical_gene_set = None


for i, path in enumerate(files, start=1):

    donor_expected = donor_from_path(path)

    print(
        f"\n[{i:02d}/{len(files):02d}] "
        f"{donor_expected}"
    )

    a = ad.read_h5ad(
        path,
        backed="r",
    )

    # ------------------------------------------------------------------
    # Basic shape
    # ------------------------------------------------------------------

    n_cells = int(a.n_obs)
    n_genes = int(a.n_vars)

    e_ids = pd.Index(
        a.obs_names.astype(str)
    )

    duplicate_envi_ids = int(
        e_ids.duplicated().sum()
    )

    # ------------------------------------------------------------------
    # Donor identity
    # ------------------------------------------------------------------

    if "BrNum" in a.obs.columns:
        donor_values = unique_nonempty(
            a.obs["BrNum"]
        )
    else:
        donor_values = []

    if "Dx" in a.obs.columns:
        dx_values = unique_nonempty(
            a.obs["Dx"]
        )
    else:
        dx_values = []

    donor_actual = (
        donor_values[0]
        if len(donor_values) == 1
        else "|".join(donor_values)
    )

    dx = (
        dx_values[0]
        if len(dx_values) == 1
        else "|".join(dx_values)
    )

    # ------------------------------------------------------------------
    # Exact annotation join
    # ------------------------------------------------------------------

    ann_d = ann_obs.loc[
        ann_obs["BrNum"] == donor_expected
    ].copy()

    a_ids = pd.Index(
        ann_d.index.astype(str)
    )

    duplicate_annotation_ids = int(
        a_ids.duplicated().sum()
    )

    overlap = e_ids.intersection(a_ids)

    join_n = len(overlap)

    join_pct = (
        100.0 * join_n / n_cells
        if n_cells > 0
        else np.nan
    )

    envi_missing_annotation = len(
        e_ids.difference(a_ids)
    )

    annotation_missing_envi = len(
        a_ids.difference(e_ids)
    )

    # ------------------------------------------------------------------
    # Coordinate audit
    # ------------------------------------------------------------------

    has_spatial = (
        "spatial" in a.obsm.keys()
    )

    spatial_finite = False
    spatial_shape = None

    if has_spatial:

        coords = np.asarray(
            a.obsm["spatial"]
        )

        spatial_shape = str(
            tuple(coords.shape)
        )

        spatial_finite = bool(
            np.isfinite(coords).all()
        )

    coordinate_max_abs_diff = np.nan

    if (
        join_n == n_cells
        and "x_centroid" in a.obs.columns
        and "y_centroid" in a.obs.columns
        and "x_centroid" in ann_d.columns
        and "y_centroid" in ann_d.columns
    ):

        ann_aligned = ann_d.loc[e_ids]

        ex = pd.to_numeric(
            a.obs["x_centroid"],
            errors="coerce",
        ).to_numpy()

        ey = pd.to_numeric(
            a.obs["y_centroid"],
            errors="coerce",
        ).to_numpy()

        ax = pd.to_numeric(
            ann_aligned["x_centroid"],
            errors="coerce",
        ).to_numpy()

        ay = pd.to_numeric(
            ann_aligned["y_centroid"],
            errors="coerce",
        ).to_numpy()

        coordinate_max_abs_diff = float(
            np.nanmax(
                np.abs(
                    np.concatenate(
                        [
                            ex - ax,
                            ey - ay,
                        ]
                    )
                )
            )
        )

    # ------------------------------------------------------------------
    # Cell-type coverage
    # ------------------------------------------------------------------

    if "cell_type_annotation" not in ann_d.columns:
        raise RuntimeError(
            "cell_type_annotation missing "
            "from annotated Xenium."
        )

    ct = normalize_string_series(
        ann_d["cell_type_annotation"]
    )

    ct_counts = ct.value_counts()

    celltype_record = {
        "BrNum": donor_expected,
        "Dx": dx,
    }

    for label in EXPECTED_CELL_TYPES:

        value = int(
            ct_counts.get(label, 0)
        )

        celltype_record[label] = value

    celltype_count_records.append(
        celltype_record
    )

    unexpected_celltypes = sorted(
        set(ct.unique())
        - set(EXPECTED_CELL_TYPES)
    )

    missing_celltypes = sorted(
        set(EXPECTED_CELL_TYPES)
        - set(ct.unique())
    )

    ambiguous_cells = int(
        ct.isin(
            [
                "Ambig/In/Endo",
                "Ambig/Oligo",
            ]
        ).sum()
    )

    ambiguous_pct = (
        100.0 * ambiguous_cells / len(ct)
        if len(ct)
        else np.nan
    )

    # ------------------------------------------------------------------
    # Layer coverage
    # ------------------------------------------------------------------

    missing_layers = []
    unexpected_layers = []

    if "layer_annotation" in ann_d.columns:

        la = normalize_string_series(
            ann_d["layer_annotation"]
        )

        layer_counts = la.value_counts()

        layer_record = {
            "BrNum": donor_expected,
            "Dx": dx,
        }

        for layer in EXPECTED_LAYERS:

            layer_record[layer] = int(
                layer_counts.get(layer, 0)
            )

        layer_count_records.append(
            layer_record
        )

        missing_layers = sorted(
            set(EXPECTED_LAYERS)
            - set(la.unique())
        )

        unexpected_layers = sorted(
            set(la.unique())
            - set(EXPECTED_LAYERS)
        )

    # ------------------------------------------------------------------
    # Gene consistency
    # ------------------------------------------------------------------

    genes = pd.Index(
        a.var_names.astype(str)
    )

    genes_unique = (
        genes.nunique() == len(genes)
    )

    if canonical_gene_names is None:

        canonical_gene_names = genes.copy()
        canonical_gene_set = set(genes)

        gene_order_identical = True
        gene_set_identical = True

    else:

        gene_order_identical = bool(
            genes.equals(
                canonical_gene_names
            )
        )

        gene_set_identical = bool(
            set(genes)
            == canonical_gene_set
        )

    # ------------------------------------------------------------------
    # Measured / imputed gene provenance
    # ------------------------------------------------------------------

    measured_n = np.nan
    imputed_n = np.nan

    if "measured_in_xenium" in a.var.columns:

        measured = (
            a.var["measured_in_xenium"]
            .astype(bool)
            .to_numpy()
        )

        measured_n = int(
            measured.sum()
        )

        imputed_n = int(
            len(measured) - measured.sum()
        )

    source_counts = {}

    if "expression_source" in a.var.columns:

        source_counts = (
            a.var["expression_source"]
            .astype(str)
            .value_counts()
            .to_dict()
        )

    # ------------------------------------------------------------------
    # Small expression integrity sample
    # ------------------------------------------------------------------

    if n_cells > 0:

        sample_n = min(
            64,
            n_cells,
        )

        row_idx = np.unique(
            np.linspace(
                0,
                n_cells - 1,
                sample_n,
                dtype=int,
            )
        )

    else:
        row_idx = np.array(
            [],
            dtype=int,
        )

    for layer in [
        "count_scale",
        "log1p",
        "native",
    ]:

        if layer not in a.layers.keys():
            continue

        try:

            # Small row sample; then spread columns in memory.
            block = np.asarray(
                a.layers[layer][
                    row_idx,
                    :
                ]
            )

            if block.ndim == 1:
                block = block.reshape(
                    1,
                    -1,
                )

            col_n = min(
                512,
                block.shape[1],
            )

            col_idx = np.unique(
                np.linspace(
                    0,
                    block.shape[1] - 1,
                    col_n,
                    dtype=int,
                )
            )

            sample = block[
                :,
                col_idx
            ]

            finite = np.isfinite(sample)

            expression_rows.append(
                {
                    "BrNum": donor_expected,
                    "Dx": dx,
                    "layer": layer,
                    "sample_cells": int(
                        sample.shape[0]
                    ),
                    "sample_genes": int(
                        sample.shape[1]
                    ),
                    "all_finite": bool(
                        finite.all()
                    ),
                    "nan_n": int(
                        np.isnan(sample).sum()
                    ),
                    "inf_n": int(
                        np.isinf(sample).sum()
                    ),
                    "min": float(
                        np.nanmin(sample)
                    ),
                    "max": float(
                        np.nanmax(sample)
                    ),
                    "mean": float(
                        np.nanmean(sample)
                    ),
                    "negative_n": int(
                        (sample < 0).sum()
                    ),
                }
            )

            del block
            del sample

        except Exception as exc:

            expression_rows.append(
                {
                    "BrNum": donor_expected,
                    "Dx": dx,
                    "layer": layer,
                    "sample_cells": np.nan,
                    "sample_genes": np.nan,
                    "all_finite": False,
                    "nan_n": np.nan,
                    "inf_n": np.nan,
                    "min": np.nan,
                    "max": np.nan,
                    "mean": np.nan,
                    "negative_n": np.nan,
                    "error": repr(exc),
                }
            )

    # ------------------------------------------------------------------
    # Donor metadata from annotated Xenium
    # ------------------------------------------------------------------

    metadata_fields = [
        "Age",
        "Sex",
        "slide_id",
        "run_date",
        "sample_id",
        "CaptureArea",
        "PNN",
        "tear",
    ]

    metadata = {}

    for col in metadata_fields:

        if col in ann_d.columns:
            metadata[col] = collapse_metadata(
                ann_d[col]
            )
        elif col in a.obs.columns:
            metadata[col] = collapse_metadata(
                a.obs[col]
            )
        else:
            metadata[col] = np.nan

    # ------------------------------------------------------------------
    # Row
    # ------------------------------------------------------------------

    row = {
        "BrNum": donor_expected,
        "Dx": dx,
        "file": str(path),
        "n_cells": n_cells,
        "n_genes": n_genes,
        "donor_in_file": donor_actual,
        "donor_matches_path": (
            donor_actual == donor_expected
        ),
        "n_dx_values": len(dx_values),
        "duplicate_envi_ids": duplicate_envi_ids,
        "annotated_cells": len(ann_d),
        "duplicate_annotation_ids":
            duplicate_annotation_ids,
        "join_n": join_n,
        "join_pct": join_pct,
        "envi_missing_annotation":
            envi_missing_annotation,
        "annotation_missing_envi":
            annotation_missing_envi,
        "has_spatial": has_spatial,
        "spatial_shape": spatial_shape,
        "spatial_finite": spatial_finite,
        "coordinate_max_abs_diff":
            coordinate_max_abs_diff,
        "n_cell_types": int(
            ct.nunique()
        ),
        "missing_cell_types":
            ";".join(missing_celltypes),
        "unexpected_cell_types":
            ";".join(unexpected_celltypes),
        "ambiguous_cells":
            ambiguous_cells,
        "ambiguous_pct":
            ambiguous_pct,
        "missing_layers":
            ";".join(missing_layers),
        "unexpected_layers":
            ";".join(unexpected_layers),
        "gene_names_unique":
            genes_unique,
        "gene_set_identical":
            gene_set_identical,
        "gene_order_identical":
            gene_order_identical,
        "measured_genes":
            measured_n,
        "imputed_genes":
            imputed_n,
        "reference_donor_overlap":
            donor_expected in reference_donors,
    }

    row.update(metadata)

    for key, value in source_counts.items():

        safe = re.sub(
            r"[^A-Za-z0-9]+",
            "_",
            key,
        ).strip("_")

        row[
            f"expression_source_{safe}"
        ] = int(value)

    donor_rows.append(row)

    a.file.close()


# ======================================================================================
# DATAFRAMES
# ======================================================================================

section("BUILD AUDIT TABLES")

donor_df = pd.DataFrame(
    donor_rows
).sort_values(
    ["Dx", "BrNum"]
)

celltype_df = pd.DataFrame(
    celltype_count_records
).sort_values(
    ["Dx", "BrNum"]
)

layer_df = pd.DataFrame(
    layer_count_records
).sort_values(
    ["Dx", "BrNum"]
)

expression_df = pd.DataFrame(
    expression_rows
)


# ======================================================================================
# GENE MANIFEST
# ======================================================================================

if files:

    first = ad.read_h5ad(
        files[0],
        backed="r",
    )

    gene_manifest = pd.DataFrame(
        {
            "gene": first.var_names.astype(str)
        }
    )

    for c in [
        "expression_source",
        "measured_in_xenium",
        "highly_variable",
        "means",
        "dispersions",
        "dispersions_norm",
    ]:

        if c in first.var.columns:
            gene_manifest[c] = (
                first.var[c]
                .to_numpy()
            )

    first.file.close()

else:
    gene_manifest = pd.DataFrame()


# ======================================================================================
# CELL-TYPE DONOR COVERAGE
# ======================================================================================

section("CELL-TYPE COVERAGE")

coverage_rows = []

if not celltype_df.empty:

    for ct in EXPECTED_CELL_TYPES:

        ntc = celltype_df[
            celltype_df["Dx"] == "NTC"
        ]

        scz = celltype_df[
            celltype_df["Dx"] == "SCZ"
        ]

        coverage_rows.append(
            {
                "cell_type": ct,
                "total_cells": int(
                    celltype_df[ct].sum()
                ),
                "donors_present": int(
                    (celltype_df[ct] > 0).sum()
                ),
                "NTC_donors_present": int(
                    (ntc[ct] > 0).sum()
                ),
                "SCZ_donors_present": int(
                    (scz[ct] > 0).sum()
                ),
                "NTC_cells": int(
                    ntc[ct].sum()
                ),
                "SCZ_cells": int(
                    scz[ct].sum()
                ),
                "min_cells_when_present": int(
                    celltype_df.loc[
                        celltype_df[ct] > 0,
                        ct,
                    ].min()
                )
                if (celltype_df[ct] > 0).any()
                else 0,
            }
        )

celltype_coverage_df = pd.DataFrame(
    coverage_rows
)

print(celltype_coverage_df.to_string(index=False))


# ======================================================================================
# LAYER DONOR COVERAGE
# ======================================================================================

section("LAYER COVERAGE")

layer_coverage_rows = []

if not layer_df.empty:

    for layer in EXPECTED_LAYERS:

        ntc = layer_df[
            layer_df["Dx"] == "NTC"
        ]

        scz = layer_df[
            layer_df["Dx"] == "SCZ"
        ]

        layer_coverage_rows.append(
            {
                "layer": layer,
                "total_cells": int(
                    layer_df[layer].sum()
                ),
                "donors_present": int(
                    (layer_df[layer] > 0).sum()
                ),
                "NTC_donors_present": int(
                    (ntc[layer] > 0).sum()
                ),
                "SCZ_donors_present": int(
                    (scz[layer] > 0).sum()
                ),
                "NTC_cells": int(
                    ntc[layer].sum()
                ),
                "SCZ_cells": int(
                    scz[layer].sum()
                ),
            }
        )

layer_coverage_df = pd.DataFrame(
    layer_coverage_rows
)

print(layer_coverage_df.to_string(index=False))


# ======================================================================================
# TECHNICAL / DEMOGRAPHIC AUDITS
# ======================================================================================

section("DIAGNOSIS AND DONOR COUNTS")

print(
    donor_df["Dx"]
    .value_counts(dropna=False)
    .to_string()
)

print("\nTotal cells by diagnosis:")

print(
    donor_df.groupby("Dx")[
        "n_cells"
    ].sum()
)


section("SEX × DIAGNOSIS")

if "Sex" in donor_df.columns:

    print(
        pd.crosstab(
            donor_df["Sex"],
            donor_df["Dx"],
            dropna=False,
        ).to_string()
    )


section("AGE BY DIAGNOSIS")

if "Age" in donor_df.columns:

    age_df = donor_df.copy()

    age_df["Age_numeric"] = pd.to_numeric(
        age_df["Age"],
        errors="coerce",
    )

    print(
        age_df.groupby("Dx")[
            "Age_numeric"
        ].agg(
            [
                "count",
                "mean",
                "std",
                "min",
                "max",
            ]
        ).to_string()
    )


# ======================================================================================
# TECHNICAL LEVELS vs DIAGNOSIS
# ======================================================================================

technical_summary = []

for col in [
    "run_date",
    "slide_id",
    "CaptureArea",
]:

    if col not in donor_df.columns:
        continue

    if donor_df[col].isna().all():
        continue

    section(f"{col} × DIAGNOSIS")

    print(
        pd.crosstab(
            donor_df[col],
            donor_df["Dx"],
            dropna=False,
        ).to_string()
    )

    ntc_levels = set(
        donor_df.loc[
            donor_df["Dx"] == "NTC",
            col,
        ]
        .dropna()
        .astype(str)
    )

    scz_levels = set(
        donor_df.loc[
            donor_df["Dx"] == "SCZ",
            col,
        ]
        .dropna()
        .astype(str)
    )

    shared = ntc_levels & scz_levels

    technical_summary.append(
        {
            "variable": col,
            "NTC_levels": len(
                ntc_levels
            ),
            "SCZ_levels": len(
                scz_levels
            ),
            "shared_levels": len(
                shared
            ),
            "shared_values":
                ";".join(sorted(shared)),
        }
    )

technical_df = pd.DataFrame(
    technical_summary
)


# ======================================================================================
# REFERENCE OVERLAP
# ======================================================================================

section("XENIUM ↔ snRNA BIOLOGICAL DONOR OVERLAP")

cohort_donors = set(
    donor_df["BrNum"]
)

overlap_donors = sorted(
    cohort_donors
    & reference_donors
)

print(
    "Cohort donors:",
    len(cohort_donors),
)

print(
    "Reference donors:",
    len(reference_donors),
)

print(
    "Overlap donors:",
    len(overlap_donors),
)

print(overlap_donors)

reference_overlap_df = pd.DataFrame(
    {
        "BrNum": sorted(
            cohort_donors
        )
    }
)

reference_overlap_df[
    "in_snrna_reference"
] = (
    reference_overlap_df[
        "BrNum"
    ].isin(
        reference_donors
    )
)


# ======================================================================================
# GLOBAL PASS / WARN / FAIL
# ======================================================================================

section("GLOBAL PRE-SPLIT STATUS")

checks = []


def add_check(name, status, detail):

    checks.append(
        {
            "check": name,
            "status": status,
            "detail": detail,
        }
    )

    print(
        f"{status:5s} | "
        f"{name:<40s} | "
        f"{detail}"
    )


# 24 donors
add_check(
    "24 ENVI donor files",
    "PASS"
    if len(donor_df) == 24
    else "FAIL",
    f"found={len(donor_df)}",
)

# 12/12 diagnosis
dx_counts = (
    donor_df["Dx"]
    .value_counts()
    .to_dict()
)

balanced_dx = (
    dx_counts.get("NTC", 0) == 12
    and dx_counts.get("SCZ", 0) == 12
)

add_check(
    "12 NTC + 12 SCZ",
    "PASS"
    if balanced_dx
    else "FAIL",
    str(dx_counts),
)

# one donor per file
ok = bool(
    donor_df[
        "donor_matches_path"
    ].all()
)

add_check(
    "Donor identity per H5AD",
    "PASS"
    if ok
    else "FAIL",
    (
        f"{int(donor_df['donor_matches_path'].sum())}"
        f"/{len(donor_df)} match"
    ),
)

# annotation join
join_ok = bool(
    np.allclose(
        donor_df["join_pct"],
        100.0,
        equal_nan=False,
    )
)

add_check(
    "ENVI ↔ annotation exact cell join",
    "PASS"
    if join_ok
    else "FAIL",
    (
        f"minimum coverage="
        f"{donor_df['join_pct'].min():.4f}%"
    ),
)

# duplicate IDs
duplicate_total = int(
    donor_df[
        "duplicate_envi_ids"
    ].sum()
    +
    donor_df[
        "duplicate_annotation_ids"
    ].sum()
)

add_check(
    "Unique cell IDs",
    "PASS"
    if duplicate_total == 0
    else "FAIL",
    f"duplicates={duplicate_total}",
)

# spatial
spatial_ok = bool(
    donor_df["has_spatial"].all()
    and donor_df[
        "spatial_finite"
    ].all()
)

add_check(
    "Spatial coordinates",
    "PASS"
    if spatial_ok
    else "FAIL",
    (
        f"{int(donor_df['spatial_finite'].sum())}"
        f"/{len(donor_df)} finite"
    ),
)

# gene set/order
gene_set_ok = bool(
    donor_df[
        "gene_set_identical"
    ].all()
)

gene_order_ok = bool(
    donor_df[
        "gene_order_identical"
    ].all()
)

add_check(
    "Same gene set across donors",
    "PASS"
    if gene_set_ok
    else "FAIL",
    (
        f"{int(donor_df['gene_set_identical'].sum())}"
        f"/{len(donor_df)}"
    ),
)

add_check(
    "Same gene order across donors",
    "PASS"
    if gene_order_ok
    else "FAIL",
    (
        f"{int(donor_df['gene_order_identical'].sum())}"
        f"/{len(donor_df)}"
    ),
)

# expected cell types
celltype_ok = bool(
    (
        donor_df[
            "unexpected_cell_types"
        ] == ""
    ).all()
)

add_check(
    "Expected Xenium cell-type vocabulary",
    "PASS"
    if celltype_ok
    else "WARN",
    (
        "12 literature-supported labels expected"
    ),
)

# expression sample finite
if not expression_df.empty:

    expression_ok = bool(
        expression_df[
            "all_finite"
        ].fillna(False).all()
    )

else:
    expression_ok = False

add_check(
    "Expression sample finite",
    "PASS"
    if expression_ok
    else "WARN",
    (
        f"{int(expression_df['all_finite'].fillna(False).sum())}"
        f"/{len(expression_df)} sampled layer checks"
        if not expression_df.empty
        else "no expression checks"
    ),
)

# reference overlap
add_check(
    "snRNA biological-donor overlap",
    "PASS"
    if len(overlap_donors) == 0
    else "WARN",
    (
        "none"
        if len(overlap_donors) == 0
        else ",".join(overlap_donors)
    ),
)


checks_df = pd.DataFrame(
    checks
)


# ======================================================================================
# WRITE OUTPUTS
# ======================================================================================

section("WRITE OUTPUTS")

donor_df.to_csv(
    OUT / "00_donor_audit.csv",
    index=False,
)

celltype_df.to_csv(
    OUT / "00_celltype_counts_by_donor.csv",
    index=False,
)

celltype_coverage_df.to_csv(
    OUT / "00_celltype_donor_coverage.csv",
    index=False,
)

layer_df.to_csv(
    OUT / "00_layer_counts_by_donor.csv",
    index=False,
)

layer_coverage_df.to_csv(
    OUT / "00_layer_donor_coverage.csv",
    index=False,
)

expression_df.to_csv(
    OUT / "00_expression_sample_qc.csv",
    index=False,
)

technical_df.to_csv(
    OUT / "00_technical_dx_overlap.csv",
    index=False,
)

reference_overlap_df.to_csv(
    OUT / "00_snrna_reference_overlap.csv",
    index=False,
)

gene_manifest.to_csv(
    OUT / "00_gene_manifest.csv",
    index=False,
)

checks_df.to_csv(
    OUT / "00_presplit_status.csv",
    index=False,
)


summary_lines = [
    "AIVC BETA PRE-SPLIT AUDIT",
    "=" * 100,
    "",
    f"ENVI donors: {len(donor_df)}",
    f"NTC: {(donor_df['Dx'] == 'NTC').sum()}",
    f"SCZ: {(donor_df['Dx'] == 'SCZ').sum()}",
    f"Total cells: {int(donor_df['n_cells'].sum())}",
    "",
    "REFERENCE OVERLAP:",
    (
        ", ".join(overlap_donors)
        if overlap_donors
        else "None"
    ),
    "",
    "CHECKS:",
]

for _, r in checks_df.iterrows():

    summary_lines.append(
        f"{r['status']:5s} | "
        f"{r['check']} | "
        f"{r['detail']}"
    )

summary_path = (
    OUT / "00_presplit_summary.txt"
)

summary_path.write_text(
    "\n".join(summary_lines)
    + "\n",
    encoding="utf-8",
)

print("\nCreated:")

for p in sorted(
    OUT.glob("00_*")
):
    print(" ", p)

print("\n")
print(summary_path.read_text())


# ======================================================================================
# CLOSE BACKED OBJECTS
# ======================================================================================

ann.file.close()
ref.file.close()

