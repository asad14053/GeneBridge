#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--canonical-folds", required=True, type=Path)
    p.add_argument("--donor", required=True)
    p.add_argument("--output-root", required=True, type=Path)

    p.add_argument("--donor-col", default="BrNum")
    p.add_argument("--spatial-key", default="spatial")

    return p.parse_args()


def make_leakage_safe(
    adata: ad.AnnData,
    spatial_key: str,
):
    # Critical: held-out genes must not survive in .raw
    adata.raw = None

    if spatial_key not in adata.obsm:
        raise RuntimeError(
            f"Missing obsm[{spatial_key!r}]"
        )

    spatial = np.asarray(
        adata.obsm[spatial_key]
    ).copy()

    # Keep physical coordinates only.
    for key in list(adata.obsm.keys()):
        del adata.obsm[key]

    adata.obsm[spatial_key] = spatial

    adata.obsp.clear()
    adata.varm.clear()
    adata.varp.clear()
    adata.uns.clear()


def main():

    args = parse_args()

    print("=" * 100)
    print("APPLY CANONICAL 3-FOLD GENE SPLIT TO N24 XENIUM")
    print("=" * 100)

    print("Input:", args.input)
    print("Canonical folds:", args.canonical_folds)
    print("Donor:", args.donor)

    if not args.input.exists():
        raise FileNotFoundError(args.input)

    if not args.canonical_folds.exists():
        raise FileNotFoundError(args.canonical_folds)

    # ------------------------------------------------------------------
    # Load canonical Br8667 gene-fold assignment
    # ------------------------------------------------------------------

    print("\n" + "=" * 100)
    print("LOAD CANONICAL FOLD ASSIGNMENT")
    print("=" * 100)

    canonical = ad.read_h5ad(
        args.canonical_folds,
        backed="r",
    )

    canonical_genes = pd.Index(
        canonical.var_names.astype(str)
    )

    if len(canonical_genes) != 300:
        raise RuntimeError(
            f"Expected 300 canonical genes; found "
            f"{len(canonical_genes)}"
        )

    if "heldout_fold" not in canonical.var.columns:
        raise RuntimeError(
            "Canonical master H5AD lacks var['heldout_fold']"
        )

    fold_assignment = pd.Series(
        canonical.var["heldout_fold"]
        .astype(int)
        .to_numpy(),
        index=canonical_genes,
        name="heldout_fold",
    )

    canonical_var = canonical.var.copy()
    canonical_var.index = canonical_genes

    print("\nCanonical fold sizes:")
    print(
        fold_assignment
        .value_counts()
        .sort_index()
        .to_string()
    )

    if set(fold_assignment.unique()) != {1, 2, 3}:
        raise RuntimeError(
            "Canonical heldout_fold must contain 1, 2, 3"
        )

    for fold in (1, 2, 3):
        if (fold_assignment == fold).sum() != 100:
            raise RuntimeError(
                f"Fold {fold} does not contain 100 genes."
            )

    # Save names for diagnostics
    for fold in (1, 2, 3):
        genes = canonical_genes[
            fold_assignment.to_numpy() == fold
        ]

        print(
            f"\nFold {fold} held-out: {len(genes)} genes"
        )

    canonical.file.close()

    # ------------------------------------------------------------------
    # Load N24 object backed
    # ------------------------------------------------------------------

    print("\n" + "=" * 100)
    print("LOAD N24")
    print("=" * 100)

    whole = ad.read_h5ad(
        args.input,
        backed="r",
    )

    if whole.n_vars != 300:
        raise RuntimeError(
            f"N24 expected 300 genes; found {whole.n_vars}"
        )

    if args.donor_col not in whole.obs.columns:
        raise RuntimeError(
            f"Missing obs[{args.donor_col!r}]"
        )

    input_genes = pd.Index(
        whole.var_names.astype(str)
    )

    if set(input_genes) != set(canonical_genes):
        missing = sorted(
            set(canonical_genes) - set(input_genes)
        )
        extra = sorted(
            set(input_genes) - set(canonical_genes)
        )

        raise RuntimeError(
            "N24 and canonical gene sets differ.\n"
            f"Missing: {missing}\n"
            f"Extra: {extra}"
        )

    donor_values = (
        whole.obs[args.donor_col]
        .astype(str)
        .str.strip()
    )

    available_donors = sorted(
        donor_values.unique()
    )

    print("N24 donors:", len(available_donors))

    if len(available_donors) != 24:
        raise RuntimeError(
            f"Expected 24 donors; found {len(available_donors)}"
        )

    if args.donor not in available_donors:
        raise RuntimeError(
            f"{args.donor} is not present in N24."
        )

    idx = np.flatnonzero(
        donor_values.to_numpy() == args.donor
    )

    print(
        f"{args.donor} cells:",
        len(idx),
    )

    if len(idx) == 0:
        raise RuntimeError(
            f"No cells for {args.donor}"
        )

    # Load ONLY this donor into RAM.
    donor = whole[idx, :].to_memory()

    donor.var_names = donor.var_names.astype(str)
    donor.obs_names = donor.obs_names.astype(str)

    # Reorder to exact canonical gene order.
    donor = donor[
        :,
        canonical_genes.tolist(),
    ].copy()

    whole.file.close()

    if donor.n_vars != 300:
        raise RuntimeError(
            f"{args.donor}: expected 300 genes."
        )

    # ------------------------------------------------------------------
    # Output directory
    # ------------------------------------------------------------------

    out = (
        args.output_root
        / args.donor
    )

    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ------------------------------------------------------------------
    # Add canonical fold information to donor.var
    # ------------------------------------------------------------------

    donor.var["heldout_fold"] = (
        fold_assignment.loc[
            donor.var_names
        ].to_numpy()
    )

    # Copy canonical gene-level balancing metrics when available.
    for column in canonical_var.columns:

        if column in donor.var.columns:
            continue

        values = canonical_var.loc[
            donor.var_names,
            column,
        ]

        donor.var[column] = values.to_numpy()

    # ------------------------------------------------------------------
    # Create three folds
    # ------------------------------------------------------------------

    manifest = {
        "donor": args.donor,
        "input": str(args.input.resolve()),
        "canonical_fold_source": str(
            args.canonical_folds.resolve()
        ),
        "n_cells": int(donor.n_obs),
        "n_genes": int(donor.n_vars),
        "folds": {},
    }

    heldout_sets = []

    for fold in (1, 2, 3):

        print("\n" + "-" * 100)
        print(f"{args.donor}: FOLD {fold}")
        print("-" * 100)

        heldout_genes = fold_assignment.index[
            fold_assignment == fold
        ].tolist()

        observed_genes = fold_assignment.index[
            fold_assignment != fold
        ].tolist()

        if len(heldout_genes) != 100:
            raise RuntimeError(
                f"Fold {fold}: expected 100 held-out genes"
            )

        if len(observed_genes) != 200:
            raise RuntimeError(
                f"Fold {fold}: expected 200 observed genes"
            )

        observed = donor[
            :,
            observed_genes,
        ].copy()

        heldout = donor[
            :,
            heldout_genes,
        ].copy()

        # Prevent hidden genes from leaking through raw/embeddings.
        make_leakage_safe(
            observed,
            args.spatial_key,
        )

        make_leakage_safe(
            heldout,
            args.spatial_key,
        )

        observed.var["split_role"] = (
            "observed_model_input"
        )

        heldout.var["split_role"] = (
            "heldout_ground_truth"
        )

        shared = {
            "donor": args.donor,
            "fold": fold,
            "canonical_fold_source": str(
                args.canonical_folds.resolve()
            ),
            "observed_gene_count": 200,
            "heldout_gene_count": 100,
            "observed_genes": observed_genes,
            "heldout_genes": heldout_genes,
        }

        observed.uns["fold_design"] = {
            **shared,
            "role": "observed_model_input",
        }

        heldout.uns["fold_design"] = {
            **shared,
            "role": "heldout_ground_truth",
            "warning": (
                "Never provide this H5AD to the "
                "imputation model."
            ),
        }

        observed_path = (
            out
            / f"fold_{fold}_observed_genes.h5ad"
        )

        heldout_path = (
            out
            / f"fold_{fold}_heldout_genes.h5ad"
        )

        observed.write_h5ad(
            observed_path,
            compression="gzip",
        )

        heldout.write_h5ad(
            heldout_path,
            compression="gzip",
        )

        # Also save explicit gene lists.
        (
            out
            / f"fold_{fold}_observed_genes.txt"
        ).write_text(
            "\n".join(observed_genes) + "\n"
        )

        (
            out
            / f"fold_{fold}_heldout_genes.txt"
        ).write_text(
            "\n".join(heldout_genes) + "\n"
        )

        manifest["folds"][str(fold)] = {
            "observed": str(
                observed_path.resolve()
            ),
            "heldout": str(
                heldout_path.resolve()
            ),
            "n_observed": 200,
            "n_heldout": 100,
        }

        heldout_sets.append(
            set(heldout_genes)
        )

        print(
            "Observed:",
            observed.shape,
        )

        print(
            "Held-out:",
            heldout.shape,
        )

    # ------------------------------------------------------------------
    # Final validation
    # ------------------------------------------------------------------

    heldout_union = set().union(
        *heldout_sets
    )

    if len(heldout_union) != 300:
        raise RuntimeError(
            "Held-out union does not equal 300 genes."
        )

    for i in range(3):
        for j in range(i + 1, 3):

            overlap = (
                heldout_sets[i]
                & heldout_sets[j]
            )

            if overlap:
                raise RuntimeError(
                    f"Folds {i+1} and {j+1} overlap: "
                    f"{sorted(overlap)}"
                )

    with (
        out / "fold_manifest.json"
    ).open("w") as handle:

        json.dump(
            manifest,
            handle,
            indent=2,
        )

    print("\n" + "=" * 100)
    print("FINAL VALIDATION")
    print("=" * 100)

    print("Donor:", args.donor)
    print("Cells:", donor.n_obs)
    print("Genes:", donor.n_vars)
    print("Fold 1 held-out:", len(heldout_sets[0]))
    print("Fold 2 held-out:", len(heldout_sets[1]))
    print("Fold 3 held-out:", len(heldout_sets[2]))
    print("Held-out union:", len(heldout_union))
    print("Pairwise overlap: 0")

    print(
        "\nFINAL STATUS: PASS — "
        "canonical 3-fold design applied to donor."
    )


if __name__ == "__main__":
    main()
