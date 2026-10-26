#!/usr/bin/env python3

"""
04_build_model_oof_celltype_pseudobulk.py

TASK 2.3 — DONOR × CELL-TYPE PSEUDOBULK
========================================

Purpose
-------
Construct donor × cell-type pseudobulk expression matrices from the
300-gene OOF-imputed predictions for all six NTC-trained models.

This adapts the LieberInstitute spatialDLPFC_SCZ_XENIUM GitHub workflow:

    code/analysis/07_cell_type_de/
        01_create_cell_type_pseudobulk_spe.R

GitHub grouping:
    var_registration = "annots"
    var_sample_id     = "BrNum"
    min_ncells        = 10

In GeneBridge, `cell_type_annotation` already contains the same 12 broad
cell-type labels as the GitHub `annots` variable.

Expression used
---------------
layers["count_scale"]

DO NOT use:
    ground_truth_counts
    ground_truth_log1p_norm

Expected cohort
---------------
24 donors
12 NTC
12 SCZ

Expected cell types
-------------------
12 cell types / donor

Expected pseudobulks
--------------------
24 donors × 12 cell types = 288 pseudobulks

Expected genes
--------------
300 OOF benchmark genes / model
"""

from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# ======================================================================================
# CONFIGURATION
# ======================================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

OUT_ROOT = (
    ROOT
    / "outputs/imputation_full/DE/task2_3_ntc_dysregulation"
    / "celltype_pseudobulk"
)

OUT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)


MODELS = {
    "VISTA": {
        "root": ROOT / "outputs/imputation_full/VISTA_N24",
        "subpath": Path("ex5_1/vista/oof_300/oof_predicted_300genes.h5ad"),
    },

    "gimVI": {
        "root": ROOT / "outputs/imputation_full/GIMVI_N24",
        "subpath": Path("ex5_1/gimvi/oof_300/oof_predicted_300genes.h5ad"),
    },

    "ENVI": {
        "root": ROOT / "outputs/imputation_full/ENVI_N24",
        "subpath": Path("ex5_1/envi/oof_300/oof_predicted_300genes.h5ad"),
    },

    "TransImp": {
        "root": ROOT / "outputs/imputation_full/TRANSIMP_N24",
        "subpath": Path("ex5_1/transimpspa/oof_300/oof_predicted_300genes.h5ad"),
    },

    "Tangram": {
        "root": ROOT / "outputs/imputation_full/TANGRAM_N24",
        "subpath": Path("ex5_1/tangram/oof_300/oof_predicted_300genes.h5ad"),
    },

    "SpaGE": {
        "root": ROOT / "outputs/imputation_full/SPAGE_N24",
        "subpath": Path("ex5_1/spage/oof_300/oof_predicted_300genes.h5ad"),
    },
}


DONORS = [
    "Br1113",
    "Br1139",
    "Br2039",
    "Br2421",
    "Br2719",
    "Br5314",
    "Br5373",
    "Br5400",
    "Br5436",
    "Br5588",
    "Br5590",
    "Br5622",
    "Br5639",
    "Br5746",
    "Br5931",
    "Br5973",
    "Br6032",
    "Br6389",
    "Br6432",
    "Br6437",
    "Br6496",
    "Br8433",
    "Br8667",
    "Br8772",
]


EXPECTED_CELL_TYPES = [
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


EXPRESSION_LAYER = "count_scale"

# Canonical Xenium technical metadata.
# Required to reproduce GitHub 07_cell_type_de covariates.
RUN_METADATA_FILE = (
    ROOT
    / "data/metadata/xenium_PI_run_metadata_24.csv"
)

MIN_CELLS = 10

EXPECTED_GENES = 300
EXPECTED_DONORS = 24
EXPECTED_CELLTYPES = 12
EXPECTED_PSEUDOBULKS = 288


# ======================================================================================
# HELPERS
# ======================================================================================

def section(title):
    print()
    print("=" * 110)
    print(title)
    print("=" * 110)


def dense_sum(X, idx):
    """
    Sum selected cells over genes.
    Works for dense or sparse matrices.
    """

    x = X[idx, :]

    if sparse.issparse(x):
        return np.asarray(
            x.sum(axis=0)
        ).ravel()

    return np.asarray(
        x
    ).sum(axis=0)


def get_scalar_metadata(obs, column, donor):
    """
    Require exactly one non-null metadata value within a donor.
    """

    if column not in obs.columns:
        return None

    values = (
        obs[column]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if len(values) == 0:
        return None

    if len(values) != 1:
        raise RuntimeError(
            f"{donor}: metadata column {column} "
            f"has multiple values: {values}"
        )

    return values[0]


# ======================================================================================
# LOAD CANONICAL SLIDE / RUN METADATA
# ======================================================================================

if not RUN_METADATA_FILE.exists():
    raise FileNotFoundError(
        f"Missing run metadata file: {RUN_METADATA_FILE}"
    )

run_meta = pd.read_csv(RUN_METADATA_FILE)

if "BrNum" not in run_meta.columns:
    raise RuntimeError(
        f"BrNum missing from {RUN_METADATA_FILE}"
    )

required_run_cols = ["BrNum", "slide_id", "run_date"]

missing_run_cols = [
    c for c in required_run_cols
    if c not in run_meta.columns
]

if missing_run_cols:
    raise RuntimeError(
        f"Run metadata missing columns: {missing_run_cols}\n"
        f"Available: {list(run_meta.columns)}"
    )

run_meta = (
    run_meta[required_run_cols]
    .drop_duplicates()
    .copy()
)

if run_meta["BrNum"].duplicated().any():
    raise RuntimeError(
        "Run metadata contains duplicated BrNum entries."
    )

RUN_META_LOOKUP = (
    run_meta
    .set_index("BrNum")
    .to_dict(orient="index")
)


# ======================================================================================
# PROCESS ONE MODEL
# ======================================================================================

def process_model(model, config):

    section(
        f"MODEL: {model}"
    )

    model_out = (
        OUT_ROOT
        / model
    )

    model_out.mkdir(
        parents=True,
        exist_ok=True
    )


    all_pb = []
    all_meta = []

    reference_genes = None

    donor_summary = []


    for donor in DONORS:

        print()
        print("-" * 90)
        print(f"{model} / {donor}")
        print("-" * 90)


        h5ad_file = (
            config["root"]
            / donor
            / config["subpath"]
        )


        if not h5ad_file.exists():
            raise FileNotFoundError(
                f"{model} {donor}: missing OOF file:\n"
                f"{h5ad_file}"
            )


        a = ad.read_h5ad(
            h5ad_file
        )


        # ==================================================================================
        # BASIC VALIDATION
        # ==================================================================================

        if a.n_vars != EXPECTED_GENES:

            raise RuntimeError(
                f"{model} {donor}: expected "
                f"{EXPECTED_GENES} genes, found {a.n_vars}"
            )


        if EXPRESSION_LAYER not in a.layers:

            raise RuntimeError(
                f"{model} {donor}: missing "
                f"layers['{EXPRESSION_LAYER}']"
            )


        required_obs = [
            "BrNum",
            "Dx",
            "Age",
            "Sex",
            "cell_type_annotation",
        ]


        missing_obs = [
            x for x in required_obs
            if x not in a.obs.columns
        ]


        if missing_obs:

            raise RuntimeError(
                f"{model} {donor}: missing obs columns: "
                f"{missing_obs}"
            )


        donor_values = (
            a.obs["BrNum"]
            .astype(str)
            .unique()
            .tolist()
        )


        if donor_values != [donor]:

            raise RuntimeError(
                f"{model} {donor}: BrNum mismatch: "
                f"{donor_values}"
            )


        # ==================================================================================
        # GENE ORDER VALIDATION
        # ==================================================================================

        genes = (
            a.var_names
            .astype(str)
            .tolist()
        )


        if reference_genes is None:

            reference_genes = genes


        elif genes != reference_genes:

            raise RuntimeError(
                f"{model} {donor}: OOF gene order "
                f"does not match first donor."
            )


        # ==================================================================================
        # CELL-TYPE VALIDATION
        # ==================================================================================

        ct = (
            a.obs["cell_type_annotation"]
            .astype(str)
        )


        observed_celltypes = set(
            ct.unique()
        )


        unexpected = (
            observed_celltypes
            - set(EXPECTED_CELL_TYPES)
        )


        if unexpected:

            raise RuntimeError(
                f"{model} {donor}: unexpected "
                f"cell-type labels: {sorted(unexpected)}"
            )


        missing_ct = (
            set(EXPECTED_CELL_TYPES)
            - observed_celltypes
        )


        if missing_ct:

            raise RuntimeError(
                f"{model} {donor}: missing "
                f"cell types: {sorted(missing_ct)}"
            )


        # ==================================================================================
        # DONOR METADATA
        # ==================================================================================

        dx = get_scalar_metadata(
            a.obs,
            "Dx",
            donor
        )

        age = get_scalar_metadata(
            a.obs,
            "Age",
            donor
        )

        sex = get_scalar_metadata(
            a.obs,
            "Sex",
            donor
        )


        # ------------------------------------------------------------------
        # Restore paper/GitHub technical metadata.
        #
        # OOF H5AD files do not carry slide_id/run_date, so these must come
        # from the canonical N24 Xenium metadata table.
        # ------------------------------------------------------------------

        if donor not in RUN_META_LOOKUP:
            raise RuntimeError(
                f"{model} {donor}: donor missing from run metadata."
            )

        slide_id = str(
            RUN_META_LOOKUP[donor]["slide_id"]
        )

        run_date = str(
            RUN_META_LOOKUP[donor]["run_date"]
        )


        # ==================================================================================
        # EXPRESSION
        # ==================================================================================

        X = a.layers[
            EXPRESSION_LAYER
        ]


        if sparse.issparse(X):

            if X.data.size > 0 and not np.isfinite(X.data).all():

                raise RuntimeError(
                    f"{model} {donor}: non-finite "
                    f"values in count_scale"
                )

        else:

            if not np.isfinite(
                np.asarray(X)
            ).all():

                raise RuntimeError(
                    f"{model} {donor}: non-finite "
                    f"values in count_scale"
                )


        # ==================================================================================
        # DONOR × CELL-TYPE PSEUDOBULK
        #
        # Equivalent conceptual grouping to GitHub:
        #
        # registration_pseudobulk(
        #     var_registration = "annots",
        #     var_sample_id = "BrNum",
        #     min_ncells = 10
        # )
        # ==================================================================================

        counts_per_ct = (
            ct
            .value_counts()
            .to_dict()
        )


        print(
            "Cells:",
            a.n_obs
        )


        print(
            "Cell types:",
            len(counts_per_ct)
        )


        print(
            "Minimum cell count:",
            min(counts_per_ct.values())
        )


        donor_pb_count = 0


        for celltype in EXPECTED_CELL_TYPES:

            mask = (
                ct.values
                == celltype
            )

            idx = np.flatnonzero(
                mask
            )

            n_cells = len(
                idx
            )


            # GitHub min_ncells = 10

            if n_cells < MIN_CELLS:

                print(
                    f"SKIP: {donor} / {celltype}: "
                    f"{n_cells} cells < {MIN_CELLS}"
                )

                continue


            pb_values = dense_sum(
                X,
                idx
            )


            pb_id = (
                f"{donor}__"
                f"{celltype.replace(' ', '_').replace('/', '_')}"
            )


            all_pb.append(
                pd.Series(
                    pb_values,
                    index=reference_genes,
                    name=pb_id,
                )
            )


            all_meta.append(
                {
                    "pseudobulk_id":
                        pb_id,

                    "BrNum":
                        donor,

                    "Dx":
                        dx,

                    "Age":
                        age,

                    "Sex":
                        sex,

                    "slide_id":
                        slide_id,

                    "run_date":
                        run_date,

                    # GitHub equivalent variable
                    "annots":
                        celltype,

                    # Keep explicit GeneBridge name too
                    "cell_type_annotation":
                        celltype,

                    "n_cells":
                        n_cells,

                    "model":
                        model,

                    "experiment":
                        "ex5_1",

                    "training":
                        "NTC",

                    "expression_source":
                        "OOF_count_scale",
                }
            )


            donor_pb_count += 1


        donor_summary.append(
            {
                "model":
                    model,

                "BrNum":
                    donor,

                "Dx":
                    dx,

                "n_cells":
                    a.n_obs,

                "cell_types_observed":
                    len(counts_per_ct),

                "pseudobulks_created":
                    donor_pb_count,

                "min_cells_celltype":
                    min(counts_per_ct.values()),

                "max_cells_celltype":
                    max(counts_per_ct.values()),
            }
        )


        del a


    # ======================================================================================
    # COMBINE MODEL
    # ======================================================================================

    pb = pd.DataFrame(
        all_pb
    )


    meta = pd.DataFrame(
        all_meta
    )


    donor_summary = pd.DataFrame(
        donor_summary
    )


    section(
        f"{model} FINAL VALIDATION"
    )


    print(
        "PB matrix:",
        pb.shape
    )

    print(
        "Metadata:",
        meta.shape
    )

    print(
        "Donors:",
        meta["BrNum"].nunique()
    )

    print(
        "Dx:"
    )

    print(
        meta[
            ["BrNum", "Dx"]
        ]
        .drop_duplicates()
        ["Dx"]
        .value_counts()
        .to_string()
    )


    print(
        "Cell types:",
        meta["annots"].nunique()
    )


    print(
        "Pseudobulks:",
        len(meta)
    )


    print(
        "Minimum cells / PB:",
        meta["n_cells"].min()
    )


    print(
        "Maximum cells / PB:",
        meta["n_cells"].max()
    )


    # ======================================================================================
    # STRICT EXPECTATIONS
    # ======================================================================================

    if pb.shape[1] != EXPECTED_GENES:

        raise RuntimeError(
            f"{model}: expected {EXPECTED_GENES} genes, "
            f"found {pb.shape[1]}"
        )


    if meta["BrNum"].nunique() != EXPECTED_DONORS:

        raise RuntimeError(
            f"{model}: expected {EXPECTED_DONORS} donors"
        )


    if meta["annots"].nunique() != EXPECTED_CELLTYPES:

        raise RuntimeError(
            f"{model}: expected {EXPECTED_CELLTYPES} cell types"
        )


    if len(meta) != EXPECTED_PSEUDOBULKS:

        raise RuntimeError(
            f"{model}: expected "
            f"{EXPECTED_PSEUDOBULKS} pseudobulks, "
            f"found {len(meta)}"
        )


    if meta["n_cells"].min() < MIN_CELLS:

        raise RuntimeError(
            f"{model}: pseudobulk below "
            f"minimum {MIN_CELLS} cells"
        )


    # ======================================================================================
    # ADD LIBRARY SIZE
    # ======================================================================================

    meta["library_size"] = (
        pb.sum(axis=1)
        .values
    )


    # ======================================================================================
    # OUTPUT PATHS
    # ======================================================================================

    pb_file = (
        model_out
        / f"{model}_NTCtrained_OOF_donor_celltype_pseudobulk_300genes.csv.gz"
    )


    meta_file = (
        model_out
        / f"{model}_NTCtrained_OOF_donor_celltype_pseudobulk_metadata.csv"
    )


    donor_file = (
        model_out
        / f"{model}_NTCtrained_OOF_donor_celltype_validation.csv"
    )


    manifest_file = (
        model_out
        / f"{model}_NTCtrained_OOF_donor_celltype_manifest.json"
    )


    # ======================================================================================
    # SAVE
    # ======================================================================================

    pb.to_csv(
        pb_file,
        compression="gzip"
    )


    meta.to_csv(
        meta_file,
        index=False
    )


    donor_summary.to_csv(
        donor_file,
        index=False
    )


    manifest = {
        "model":
            model,

        "experiment":
            "ex5_1",

        "training":
            "NTC",

        "expression_layer":
            EXPRESSION_LAYER,

        "grouping":
            "BrNum x cell_type_annotation",

        "github_equivalent_registration_variable":
            "annots",

        "min_ncells":
            MIN_CELLS,

        "donors":
            int(meta["BrNum"].nunique()),

        "cell_types":
            int(meta["annots"].nunique()),

        "pseudobulks":
            int(len(meta)),

        "genes":
            int(pb.shape[1]),

        "cell_type_labels":
            EXPECTED_CELL_TYPES,
    }


    with open(
        manifest_file,
        "w"
    ) as handle:

        json.dump(
            manifest,
            handle,
            indent=2
        )


    print()
    print(
        "PB:",
        pb_file
    )

    print(
        "Metadata:",
        meta_file
    )

    print(
        "Validation:",
        donor_file
    )

    print(
        "Manifest:",
        manifest_file
    )


    print()
    print(
        f"STATUS: {model} PASS"
    )


# ======================================================================================
# MAIN
# ======================================================================================

section(
    "TASK 2.3 — BUILD DONOR × CELL-TYPE OOF PSEUDOBULK"
)


for model, config in MODELS.items():

    process_model(
        model,
        config
    )


section(
    "FINAL STATUS"
)

print(
    "PASS — ALL SIX MODELS COMPLETE"
)

print(
    f"Expected/model: "
    f"{EXPECTED_DONORS} donors × "
    f"{EXPECTED_CELLTYPES} cell types = "
    f"{EXPECTED_PSEUDOBULKS} pseudobulks"
)

print(
    f"Genes/model: {EXPECTED_GENES}"
)

print(
    f"Expression layer: {EXPRESSION_LAYER}"
)
