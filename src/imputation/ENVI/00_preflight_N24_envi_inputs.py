#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

FOLD_ROOT = (
    ROOT
    / "data/processed/imputation_full/gene_folds_200_100"
)

REF_ROOT = (
    ROOT
    / "data/processed/imputation_full/snrnaseq_target_references"
)

FULL_REF = (
    ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated/"
      "huuki_snrna_reference_full_allgenes.h5ad"
)

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


REFERENCES = {
    "ex5_1_full10": FULL_REF,

    "ex5_Br6432":
        REF_ROOT
        / "Br6432"
        / "seq_data_huuki_snrna_Br6432.h5ad",

    "ex5_3_Br6432":
        REF_ROOT
        / "Br6432"
        / "seq_data_huuki_snrna_9brains_excluding_Br6432.h5ad",

    "ex5_Br8667":
        REF_ROOT
        / "Br8667"
        / "seq_data_huuki_snrna_Br8667.h5ad",

    "ex5_3_Br8667":
        REF_ROOT
        / "Br8667"
        / "seq_data_huuki_snrna_9brains_excluding_Br8667.h5ad",
}


def row_sums(x):
    if sparse.issparse(x):
        return np.asarray(
            x.sum(axis=1)
        ).ravel()

    return np.asarray(x).sum(axis=1)


def read_gene_set(path):
    a = ad.read_h5ad(
        path,
        backed="r",
    )

    genes = pd.Index(
        a.var_names.astype(str)
    )

    if not genes.is_unique:
        a.file.close()
        raise RuntimeError(
            f"Duplicate genes: {path}"
        )

    result = set(genes)

    a.file.close()

    return result


# ============================================================
# REFERENCES
# ============================================================

print("=" * 100)
print("1. REFERENCE AUDIT")
print("=" * 100)

reference_gene_sets = {}

for name, path in REFERENCES.items():

    if not path.exists():
        raise FileNotFoundError(
            f"{name}: {path}"
        )

    a = ad.read_h5ad(
        path,
        backed="r",
    )

    genes = pd.Index(
        a.var_names.astype(str)
    )

    if not genes.is_unique:
        raise RuntimeError(
            f"{name}: duplicate reference genes."
        )

    brains = []

    if "BrNum" in a.obs.columns:
        brains = sorted(
            a.obs["BrNum"]
            .astype(str)
            .unique()
            .tolist()
        )

    reference_gene_sets[name] = set(
        genes
    )

    print()
    print(name)
    print("-" * 80)
    print("Path  :", path)
    print("Shape :", a.shape)
    print("Genes :", len(genes))

    if brains:
        print("Brains:", brains)
        print("N brain:", len(brains))

    # Cohort-specific assertions
    if name == "ex5_1_full10":
        if len(brains) != 10:
            raise RuntimeError(
                "Ex-5.1 must contain 10 brains."
            )

    if name == "ex5_Br6432":
        if brains != ["Br6432"]:
            raise RuntimeError(
                "Br6432 Ex-5 reference is incorrect."
            )

    if name == "ex5_Br8667":
        if brains != ["Br8667"]:
            raise RuntimeError(
                "Br8667 Ex-5 reference is incorrect."
            )

    if name == "ex5_3_Br6432":
        if len(brains) != 9:
            raise RuntimeError(
                "Br6432 Ex-5.3 must contain 9 brains."
            )
        if "Br6432" in brains:
            raise RuntimeError(
                "Br6432 remains in its Ex-5.3 reference."
            )

    if name == "ex5_3_Br8667":
        if len(brains) != 9:
            raise RuntimeError(
                "Br8667 Ex-5.3 must contain 9 brains."
            )
        if "Br8667" in brains:
            raise RuntimeError(
                "Br8667 remains in its Ex-5.3 reference."
            )

    a.file.close()


# ============================================================
# ALL REFERENCE GENE UNIVERSES
# ============================================================

print()
print("=" * 100)
print("2. REFERENCE GENE-UNIVERSE CONSISTENCY")
print("=" * 100)

reference_names = list(
    reference_gene_sets
)

base_name = reference_names[0]
base_genes = reference_gene_sets[
    base_name
]

for name in reference_names[1:]:

    current = reference_gene_sets[
        name
    ]

    same = (
        current
        == base_genes
    )

    print(
        f"{base_name} vs {name}:",
        same,
    )

    if not same:
        print(
            "Only base:",
            len(
                base_genes
                - current
            )
        )

        print(
            "Only current:",
            len(
                current
                - base_genes
            )
        )

        raise RuntimeError(
            "Reference gene universes differ."
        )


print()
print(
    "REFERENCE GENE UNIVERSE: PASS"
)


# ============================================================
# N24 FOLD AUDIT
# ============================================================

print()
print("=" * 100)
print("3. N24 DONOR × FOLD AUDIT")
print("=" * 100)

canonical_fold_heldout = {}

rows = []

for donor in DONORS:

    donor_dir = (
        FOLD_ROOT
        / donor
    )

    if not donor_dir.exists():
        raise FileNotFoundError(
            donor_dir
        )

    donor_heldout_sets = []

    canonical_cells = None

    for fold in (1, 2, 3):

        observed_path = (
            donor_dir
            / f"fold_{fold}_observed_genes.h5ad"
        )

        heldout_path = (
            donor_dir
            / f"fold_{fold}_heldout_genes.h5ad"
        )

        if not observed_path.exists():
            raise FileNotFoundError(
                observed_path
            )

        if not heldout_path.exists():
            raise FileNotFoundError(
                heldout_path
            )

        observed = ad.read_h5ad(
            observed_path
        )

        heldout = ad.read_h5ad(
            heldout_path
        )

        observed.obs_names = (
            observed.obs_names.astype(str)
        )

        heldout.obs_names = (
            heldout.obs_names.astype(str)
        )

        observed_genes = set(
            observed.var_names.astype(str)
        )

        heldout_genes = set(
            heldout.var_names.astype(str)
        )

        # --------------------------------------------
        # Dimensions
        # --------------------------------------------

        if observed.n_vars != 200:
            raise RuntimeError(
                f"{donor} Fold {fold}: "
                f"{observed.n_vars} observed genes."
            )

        if heldout.n_vars != 100:
            raise RuntimeError(
                f"{donor} Fold {fold}: "
                f"{heldout.n_vars} heldout genes."
            )

        # --------------------------------------------
        # No overlap
        # --------------------------------------------

        overlap = (
            observed_genes
            & heldout_genes
        )

        if overlap:
            raise RuntimeError(
                f"{donor} Fold {fold}: "
                f"{len(overlap)} observed/heldout "
                "gene overlap."
            )

        # --------------------------------------------
        # 300 total
        # --------------------------------------------

        if len(
            observed_genes
            | heldout_genes
        ) != 300:

            raise RuntimeError(
                f"{donor} Fold {fold}: "
                "observed + heldout != 300 genes."
            )

        # --------------------------------------------
        # Cell universe
        # --------------------------------------------

        observed_cells = (
            observed.obs_names.to_numpy()
        )

        heldout_cells = (
            heldout.obs_names.to_numpy()
        )

        if not np.array_equal(
            observed_cells,
            heldout_cells,
        ):
            raise RuntimeError(
                f"{donor} Fold {fold}: "
                "observed/heldout cells differ."
            )

        if canonical_cells is None:
            canonical_cells = (
                observed_cells.copy()
            )

        elif not np.array_equal(
            canonical_cells,
            observed_cells,
        ):
            raise RuntimeError(
                f"{donor}: cell universe differs "
                "across folds."
            )

        # --------------------------------------------
        # Blank-cell check
        # TARGET INPUT ONLY
        # --------------------------------------------

        totals = row_sums(
            observed.X
        )

        blank = int(
            np.sum(
                totals <= 0
            )
        )

        nonfinite = int(
            np.sum(
                ~np.isfinite(
                    totals
                )
            )
        )

        if blank != 0:
            raise RuntimeError(
                f"{donor} Fold {fold}: "
                f"{blank} blank observed cells remain."
            )

        if nonfinite != 0:
            raise RuntimeError(
                f"{donor} Fold {fold}: "
                f"{nonfinite} nonfinite libraries."
            )

        # --------------------------------------------
        # Reference coverage
        # --------------------------------------------

        target_300 = (
            observed_genes
            | heldout_genes
        )

        coverage = {}

        for ref_name, ref_genes in (
            reference_gene_sets.items()
        ):

            missing = (
                target_300
                - ref_genes
            )

            coverage[ref_name] = (
                len(missing)
            )

            if missing:
                raise RuntimeError(
                    f"{donor} Fold {fold}: "
                    f"{len(missing)} genes missing "
                    f"from {ref_name}: "
                    f"{sorted(missing)[:20]}"
                )

        donor_heldout_sets.append(
            heldout_genes
        )

        # --------------------------------------------
        # Canonical fold identities across donors
        # --------------------------------------------

        if fold not in canonical_fold_heldout:

            canonical_fold_heldout[
                fold
            ] = heldout_genes.copy()

        elif (
            canonical_fold_heldout[
                fold
            ]
            != heldout_genes
        ):

            raise RuntimeError(
                f"{donor} Fold {fold}: "
                "heldout gene set differs "
                "from canonical donor fold."
            )

        rows.append({
            "donor":
                donor,

            "fold":
                fold,

            "cells":
                observed.n_obs,

            "observed_genes":
                observed.n_vars,

            "heldout_genes":
                heldout.n_vars,

            "blank_observed_cells":
                blank,

            "missing_reference_genes":
                max(
                    coverage.values()
                ),
        })

        observed.file.close() if observed.isbacked else None
        heldout.file.close() if heldout.isbacked else None

        del observed
        del heldout

    # --------------------------------------------
    # Across 3 folds
    # --------------------------------------------

    union = set().union(
        *donor_heldout_sets
    )

    if len(union) != 300:
        raise RuntimeError(
            f"{donor}: heldout union is "
            f"{len(union)}, expected 300."
        )

    for i in range(3):

        for j in range(
            i + 1,
            3,
        ):

            overlap = (
                donor_heldout_sets[i]
                & donor_heldout_sets[j]
            )

            if overlap:
                raise RuntimeError(
                    f"{donor}: heldout folds "
                    f"{i+1}/{j+1} overlap."
                )

    print(
        f"{donor}: "
        f"cells={len(canonical_cells):,} | "
        "3 folds PASS | "
        "300/300 reference genes PASS"
    )


# ============================================================
# FINAL SUMMARY
# ============================================================

df = pd.DataFrame(
    rows
)

print()
print("=" * 100)
print("4. SUMMARY")
print("=" * 100)

print(
    "Donors:",
    df["donor"].nunique(),
)

print(
    "Donor-fold target inputs:",
    len(df),
)

print(
    "Observed genes per run:",
    sorted(
        df[
            "observed_genes"
        ].unique()
    ),
)

print(
    "Heldout genes per run:",
    sorted(
        df[
            "heldout_genes"
        ].unique()
    ),
)

print(
    "Remaining blank observed cells:",
    int(
        df[
            "blank_observed_cells"
        ].sum()
    ),
)

print(
    "Maximum missing reference genes:",
    int(
        df[
            "missing_reference_genes"
        ].max()
    ),
)

print()
print("Expected primary Ex-5.1 runs:")
print(
    f"  {df['donor'].nunique()} donors × 3 folds "
    f"= {len(df)} ENVI runs"
)

print()
print("Matched-reference experiments:")
print("  Br6432 Ex-5   : READY")
print("  Br6432 Ex-5.3 : READY")
print("  Br8667 Ex-5   : READY")
print("  Br8667 Ex-5.3 : READY")

print()
print("=" * 100)
print("FINAL STATUS: PASS — ALL N24 ENVI INPUTS ARE READY")
print("=" * 100)
