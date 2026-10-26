#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


PROJECT_ROOT_DEFAULT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)


def dense32(x):
    if sp.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def sparse32(x):
    if sp.issparse(x):
        return x.astype(np.float32).tocsr()
    return sp.csr_matrix(np.asarray(x, dtype=np.float32))


def check_finite(name, x):
    if sp.issparse(x):
        ok = np.isfinite(x.data).all()
    else:
        ok = np.isfinite(x).all()

    if not ok:
        raise ValueError(f"{name} contains NaN/Inf")


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--project-root",
        default=str(PROJECT_ROOT_DEFAULT),
    )

    parser.add_argument(
        "--donor",
        required=True,
    )

    parser.add_argument(
        "--experiment",
        required=True,
    )

    args = parser.parse_args()

    root = Path(args.project_root)
    donor = args.donor
    experiment = args.experiment

    # ------------------------------------------------------------------
    # Canonical 300-gene Xenium ordering
    # ------------------------------------------------------------------

    canonical_h5ad = (
        root
        / "data/processed/xenium"
        / "xenium_N24_layer_celltype_annotated.h5ad"
    )

    if not canonical_h5ad.exists():
        raise FileNotFoundError(
            f"Canonical Xenium file missing: {canonical_h5ad}"
        )

    canonical = ad.read_h5ad(
        canonical_h5ad,
        backed="r",
    )

    canonical_genes = list(
        map(str, canonical.var_names)
    )

    canonical.file.close()

    if len(canonical_genes) != 300:
        raise ValueError(
            f"Canonical Xenium has {len(canonical_genes)} genes, "
            "expected 300."
        )

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------

    fold_root = (
        root
        / "data/processed/imputation_full"
        / "gene_folds_200_100"
        / donor
    )

    transimp_root = (
        root
        / "outputs/imputation_full/TRANSIMP_N24"
        / donor
        / experiment
        / "transimpspa"
    )

    out_dir = transimp_root / "oof_300"
    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_h5ad = (
        out_dir
        / "oof_predicted_300genes.h5ad"
    )

    summary_json = (
        out_dir
        / "oof_stitch_summary.json"
    )

    print("=" * 110)
    print("STITCH TransImpSpa 3-FOLD OOF PREDICTIONS")
    print("=" * 110)
    print("Donor      :", donor)
    print("Experiment :", experiment)
    print("Fold root  :", fold_root)
    print("TransImpSpa root  :", transimp_root)
    print("Output     :", output_h5ad)
    print("=" * 110)

    # ------------------------------------------------------------------
    # Containers
    # ------------------------------------------------------------------

    gene_order_foldwise = []

    pred_count_parts = []
    pred_log1p_parts = []
    pred_native_parts = []

    truth_counts_parts = []
    truth_log_parts = []

    var_parts = []

    reference_cells = None
    reference_obs = None
    reference_spatial = None

    scale_factor_by_fold = {}

    fold_summary = []

    # ------------------------------------------------------------------
    # Read folds
    # ------------------------------------------------------------------

    for fold in (1, 2, 3):

        print()
        print("=" * 110)
        print(f"FOLD {fold}")
        print("=" * 110)

        fold_output = (
            transimp_root
            / f"fold_{fold}"
        )

        pred_path = (
            fold_output
            / "predicted_heldout_genes.h5ad"
        )

        truth_path = (
            fold_root
            / f"fold_{fold}_heldout_genes.h5ad"
        )

        complete_flag = (
            fold_output
            / "complete.flag"
        )

        run_complete_flag = (
            fold_output
            / "run_complete.flag"
        )

        for path in [
            pred_path,
            truth_path,
            complete_flag,
            run_complete_flag,
        ]:
            if not path.exists():
                raise FileNotFoundError(
                    f"Missing required input: {path}"
                )

        print("Prediction:", pred_path)
        print("Truth     :", truth_path)

        pred = ad.read_h5ad(pred_path)
        truth = ad.read_h5ad(truth_path)

        print("Prediction shape:", pred.shape)
        print("Truth shape     :", truth.shape)

        # --------------------------------------------------------------
        # Structural validation
        # --------------------------------------------------------------

        if pred.shape != truth.shape:
            raise ValueError(
                f"Prediction/truth shape mismatch: "
                f"{pred.shape} vs {truth.shape}"
            )

        if pred.n_vars != 100:
            raise ValueError(
                f"Fold {fold}: expected 100 genes, "
                f"found {pred.n_vars}"
            )

        if not pred.obs_names.equals(
            truth.obs_names
        ):
            raise ValueError(
                f"Fold {fold}: cell ordering mismatch"
            )

        pred_genes = set(
            map(str, pred.var_names)
        )

        truth_genes = set(
            map(str, truth.var_names)
        )

        if pred_genes != truth_genes:
            raise ValueError(
                f"Fold {fold}: prediction/truth gene sets differ"
            )

        # Ground-truth ordering is authoritative within fold.
        genes = list(
            map(str, truth.var_names)
        )

        pred = pred[:, genes].copy()
        truth = truth[:, genes].copy()

        # --------------------------------------------------------------
        # Required layers
        # --------------------------------------------------------------

        for layer in [
            "count_scale",
            "log1p",
            "native",
        ]:
            if layer not in pred.layers:
                raise ValueError(
                    f"Fold {fold}: prediction missing layer "
                    f"{layer!r}"
                )

        for layer in [
            "counts",
            "log1p_norm",
        ]:
            if layer not in truth.layers:
                raise ValueError(
                    f"Fold {fold}: truth missing layer "
                    f"{layer!r}"
                )

        # --------------------------------------------------------------
        # Cell universe across folds
        # --------------------------------------------------------------

        if reference_cells is None:

            reference_cells = (
                truth.obs_names.copy()
            )

            # Use Xenium truth metadata as canonical obs.
            reference_obs = truth.obs.copy()

            if "spatial" in truth.obsm:
                reference_spatial = np.asarray(
                    truth.obsm["spatial"]
                ).copy()

        else:

            if not truth.obs_names.equals(
                reference_cells
            ):
                raise ValueError(
                    f"Fold {fold}: cells/order differ "
                    "across folds"
                )

        # --------------------------------------------------------------
        # Preserve fold-specific TransImpSpa scale factor
        # --------------------------------------------------------------

        # TransImpSpa converts its native prediction back to count
        # scale directly from the observed Xenium library size.
        # No separate per-cell calibration-factor column is required.

        # --------------------------------------------------------------
        # Pull prediction representations
        # --------------------------------------------------------------

        pred_count = dense32(
            pred.layers["count_scale"]
        )

        pred_log1p = dense32(
            pred.layers["log1p"]
        )

        pred_native = dense32(
            pred.layers["native"]
        )

        # --------------------------------------------------------------
        # Pull measured ground truth
        # --------------------------------------------------------------

        truth_counts = sparse32(
            truth.layers["counts"]
        )

        truth_log = sparse32(
            truth.layers["log1p_norm"]
        )

        check_finite(
            f"{donor}/{experiment}/fold{fold}/pred_count",
            pred_count,
        )

        check_finite(
            f"{donor}/{experiment}/fold{fold}/pred_log1p",
            pred_log1p,
        )

        check_finite(
            f"{donor}/{experiment}/fold{fold}/pred_native",
            pred_native,
        )

        check_finite(
            f"{donor}/{experiment}/fold{fold}/truth_counts",
            truth_counts,
        )

        check_finite(
            f"{donor}/{experiment}/fold{fold}/truth_log",
            truth_log,
        )

        pred_count_parts.append(
            pred_count
        )

        pred_log1p_parts.append(
            pred_log1p
        )

        pred_native_parts.append(
            pred_native
        )

        truth_counts_parts.append(
            truth_counts
        )

        truth_log_parts.append(
            truth_log
        )

        gene_order_foldwise.extend(
            genes
        )

        fold_var = truth.var.copy()

        fold_var["oof_fold"] = fold
        fold_var["oof_role"] = (
            "heldout_prediction"
        )

        var_parts.append(
            fold_var
        )

        fold_summary.append(
            {
                "fold": fold,
                "n_cells": int(pred.n_obs),
                "n_genes": int(pred.n_vars),
                "prediction_count_scale_min": float(
                    pred_count.min()
                ),
                "prediction_count_scale_max": float(
                    pred_count.max()
                ),
                "prediction_count_scale_mean": float(
                    pred_count.mean()
                ),
                "truth_counts_sum": float(
                    truth_counts.sum()
                ),
            }
        )

        print(
            f"Fold {fold}: PASS"
        )

    # ------------------------------------------------------------------
    # Validate fold union
    # ------------------------------------------------------------------

    if len(gene_order_foldwise) != 300:
        raise ValueError(
            f"Stitched fold list has "
            f"{len(gene_order_foldwise)} genes"
        )

    if len(
        set(gene_order_foldwise)
    ) != 300:
        raise ValueError(
            "Genes overlap between held-out folds"
        )

    if set(gene_order_foldwise) != set(
        canonical_genes
    ):
        missing = sorted(
            set(canonical_genes)
            - set(gene_order_foldwise)
        )

        extra = sorted(
            set(gene_order_foldwise)
            - set(canonical_genes)
        )

        raise ValueError(
            "Fold union != canonical 300 genes.\n"
            f"Missing: {missing}\n"
            f"Extra: {extra}"
        )

    # ------------------------------------------------------------------
    # Concatenate in fold order
    # ------------------------------------------------------------------

    pred_count = np.concatenate(
        pred_count_parts,
        axis=1,
    )

    pred_log1p = np.concatenate(
        pred_log1p_parts,
        axis=1,
    )

    pred_native = np.concatenate(
        pred_native_parts,
        axis=1,
    )

    truth_counts = sp.hstack(
        truth_counts_parts,
        format="csr",
    )

    truth_log = sp.hstack(
        truth_log_parts,
        format="csr",
    )

    var = pd.concat(
        var_parts,
        axis=0,
    )

    # ------------------------------------------------------------------
    # Reorder into canonical Xenium 300-gene ordering
    # ------------------------------------------------------------------

    position = {
        gene: idx
        for idx, gene
        in enumerate(gene_order_foldwise)
    }

    reorder = np.asarray(
        [
            position[g]
            for g in canonical_genes
        ],
        dtype=int,
    )

    pred_count = (
        pred_count[:, reorder]
    )

    pred_log1p = (
        pred_log1p[:, reorder]
    )

    pred_native = (
        pred_native[:, reorder]
    )

    truth_counts = (
        truth_counts[:, reorder]
        .tocsr()
    )

    truth_log = (
        truth_log[:, reorder]
        .tocsr()
    )

    var = (
        var.loc[canonical_genes]
        .copy()
    )

    # ------------------------------------------------------------------
    # Add fold-specific scale factors
    # ------------------------------------------------------------------

    for fold, values in (
        scale_factor_by_fold.items()
    ):

        reference_obs[
            f"envi_count_scale_factor_fold{fold}"
        ] = values

    # ------------------------------------------------------------------
    # Build stitched AnnData
    #
    # Primary .X = count-space OOF prediction
    # ------------------------------------------------------------------

    out = ad.AnnData(
        X=pred_count.astype(
            np.float32,
            copy=False,
        ),
        obs=reference_obs,
        var=var,
    )

    # Explicitly preserve all important representations.

    out.layers["count_scale"] = (
        pred_count.astype(
            np.float32,
            copy=False,
        )
    )

    out.layers["log1p"] = (
        pred_log1p.astype(
            np.float32,
            copy=False,
        )
    )

    out.layers["native"] = (
        pred_native.astype(
            np.float32,
            copy=False,
        )
    )

    out.layers[
        "ground_truth_counts"
    ] = truth_counts

    out.layers[
        "ground_truth_log1p_norm"
    ] = truth_log

    if reference_spatial is not None:
        out.obsm["spatial"] = (
            reference_spatial
        )

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    out.uns["oof_stitch"] = {
        "donor": donor,
        "experiment": experiment,
        "n_folds": 3,
        "genes_per_fold": 100,
        "total_genes": 300,
        "primary_X": (
            "TransImpSpa held-out prediction, "
            "count_scale layer"
        ),
        "prediction_layers": [
            "count_scale",
            "log1p",
            "native",
        ],
        "truth_layers": [
            "ground_truth_counts",
            "ground_truth_log1p_norm",
        ],
        "gene_order": (
            "canonical xenium_N24_layer_celltype_annotated.h5ad"
        ),
        "cell_order": (
            "held-out Xenium fold truth"
        ),
    }

    # ------------------------------------------------------------------
    # Final validation
    # ------------------------------------------------------------------

    if out.n_vars != 300:
        raise ValueError(
            f"Final n_vars={out.n_vars}, expected 300"
        )

    if list(out.var_names) != canonical_genes:
        raise ValueError(
            "Final gene order is not canonical"
        )

    fold_counts = (
        out.var["oof_fold"]
        .value_counts()
        .sort_index()
        .to_dict()
    )

    if fold_counts != {
        1: 100,
        2: 100,
        3: 100,
    }:
        raise ValueError(
            f"Incorrect fold counts: {fold_counts}"
        )

    check_finite(
        "final count_scale",
        out.X,
    )

    check_finite(
        "final ground_truth_counts",
        out.layers[
            "ground_truth_counts"
        ],
    )

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------

    print()
    print("=" * 110)
    print("WRITING")
    print("=" * 110)

    out.write_h5ad(
        output_h5ad,
        compression="gzip",
    )

    summary = {
        "donor": donor,
        "experiment": experiment,
        "n_cells": int(out.n_obs),
        "n_genes": int(out.n_vars),
        "fold_counts": {
            str(k): int(v)
            for k, v in fold_counts.items()
        },
        "output": str(output_h5ad),
        "folds": fold_summary,
        "status": "PASS",
    }

    with open(
        summary_json,
        "w",
    ) as handle:
        json.dump(
            summary,
            handle,
            indent=2,
        )

    print("Shape                    :", out.shape)
    print("Primary X                : TransImpSpa count_scale")
    print("Prediction layers        :", [
        "count_scale",
        "log1p",
        "native",
    ])
    print("Truth layers             :", [
        "ground_truth_counts",
        "ground_truth_log1p_norm",
    ])
    print("Fold counts              :", fold_counts)
    print("Output                   :", output_h5ad)
    print("Summary                  :", summary_json)

    print()
    print(
        "FINAL STATUS: PASS"
    )


if __name__ == "__main__":
    main()
