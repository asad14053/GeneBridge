#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

DEFAULT_FOLD_ROOT = (
    ROOT
    / "data/processed/imputation_full/gene_folds_200_100"
)


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--donor",
        required=True,
        help="Donor ID, e.g. Br8667",
    )

    p.add_argument(
        "--fold-root",
        type=Path,
        default=DEFAULT_FOLD_ROOT,
    )

    return p.parse_args()


def row_sums(matrix) -> np.ndarray:
    """Return one expression-library total per cell."""

    if sparse.issparse(matrix):
        return np.asarray(
            matrix.sum(axis=1)
        ).ravel()

    return np.asarray(
        matrix
    ).sum(axis=1)


def matrix_values(matrix) -> np.ndarray:

    if sparse.issparse(matrix):
        return np.asarray(matrix.data)

    return np.asarray(matrix).ravel()


def validate_expression(
    name: str,
    matrix,
) -> None:

    values = matrix_values(matrix)

    if values.size == 0:
        return

    nonfinite = int(
        np.sum(
            ~np.isfinite(values)
        )
    )

    negative = int(
        np.sum(
            values < 0
        )
    )

    if nonfinite:
        raise ValueError(
            f"{name}: found {nonfinite} "
            "nonfinite expression values."
        )

    if negative:
        raise ValueError(
            f"{name}: found {negative} "
            "negative expression values."
        )


def write_atomic(
    adata: ad.AnnData,
    destination: Path,
) -> None:

    temporary = destination.with_name(
        destination.stem
        + ".tmp.h5ad"
    )

    if temporary.exists():
        temporary.unlink()

    adata.write_h5ad(
        temporary,
        compression="gzip",
    )

    os.replace(
        temporary,
        destination,
    )


def main():

    args = parse_args()

    donor = args.donor

    fold_dir = (
        args.fold_root
        / donor
    )

    if not fold_dir.exists():
        raise FileNotFoundError(
            fold_dir
        )

    backup_dir = (
        fold_dir
        / "backup_before_common_nonempty_filter"
    )

    common_mask_path = (
        fold_dir
        / "common_nonempty_cells.csv"
    )

    removed_path = (
        fold_dir
        / "removed_blank_cells.csv"
    )

    audit_path = (
        fold_dir
        / "blank_cell_audit_all_cells.csv"
    )

    manifest_path = (
        fold_dir
        / "blank_cell_filter_manifest.json"
    )

    def active_path(
        fold: int,
        kind: str,
    ) -> Path:

        return (
            fold_dir
            / f"fold_{fold}_{kind}_genes.h5ad"
        )

    def backup_path(
        fold: int,
        kind: str,
    ) -> Path:

        return (
            backup_dir
            / f"fold_{fold}_{kind}_genes.h5ad"
        )

    # ============================================================
    # Header
    # ============================================================

    print("=" * 100)
    print("COMMON NONEMPTY-CELL FILTER")
    print("=" * 100)

    print("Donor:", donor)
    print("Fold directory:", fold_dir)

    # ============================================================
    # Back up original fold files
    # ============================================================

    print("\n" + "=" * 100)
    print("CREATE / VERIFY BACKUPS")
    print("=" * 100)

    backup_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for fold in (1, 2, 3):

        for kind in (
            "observed",
            "heldout",
        ):

            source = active_path(
                fold,
                kind,
            )

            backup = backup_path(
                fold,
                kind,
            )

            if not source.exists():
                raise FileNotFoundError(
                    f"Missing fold file: {source}"
                )

            if backup.exists():

                print(
                    "Backup already exists:",
                    backup.name,
                )

            else:

                print(
                    "Creating backup:",
                    backup.name,
                )

                shutil.copy2(
                    source,
                    backup,
                )

    # ============================================================
    # Load original backups and construct common mask
    # ============================================================

    observed_data = {}
    heldout_data = {}
    observed_totals = {}

    canonical_cells = None
    common_keep = None

    per_fold_summary = {}

    print("\n" + "=" * 100)
    print("AUDIT ALL THREE FOLDS")
    print("=" * 100)

    for fold in (1, 2, 3):

        observed = ad.read_h5ad(
            backup_path(
                fold,
                "observed",
            )
        )

        heldout = ad.read_h5ad(
            backup_path(
                fold,
                "heldout",
            )
        )

        observed.obs_names = (
            observed.obs_names.astype(str)
        )

        heldout.obs_names = (
            heldout.obs_names.astype(str)
        )

        # --------------------------------------------------------
        # Basic dimensions
        # --------------------------------------------------------

        if observed.n_vars != 200:

            raise ValueError(
                f"{donor} Fold {fold}: "
                f"expected 200 observed genes, "
                f"found {observed.n_vars}."
            )

        if heldout.n_vars != 100:

            raise ValueError(
                f"{donor} Fold {fold}: "
                f"expected 100 held-out genes, "
                f"found {heldout.n_vars}."
            )

        # --------------------------------------------------------
        # Cell IDs must match exactly
        # --------------------------------------------------------

        if not np.array_equal(
            observed.obs_names.to_numpy(),
            heldout.obs_names.to_numpy(),
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "observed and held-out "
                "cell order differs."
            )

        validate_expression(
            f"{donor} Fold {fold} observed",
            observed.X,
        )

        validate_expression(
            f"{donor} Fold {fold} heldout",
            heldout.X,
        )

        cells = (
            observed.obs_names
            .to_numpy()
        )

        # IMPORTANT:
        # Cell filtering is based ONLY on the 200 observed genes.
        totals = row_sums(
            observed.X
        )

        if canonical_cells is None:

            canonical_cells = (
                cells.copy()
            )

            common_keep = np.ones(
                observed.n_obs,
                dtype=bool,
            )

        elif not np.array_equal(
            canonical_cells,
            cells,
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "cell IDs/order differ "
                "from Fold 1."
            )

        finite = np.isfinite(
            totals
        )

        nonempty = (
            totals > 0
        )

        fold_keep = (
            finite
            &
            nonempty
        )

        # A cell must be usable in ALL three folds.
        common_keep &= fold_keep

        observed_data[fold] = (
            observed
        )

        heldout_data[fold] = (
            heldout
        )

        observed_totals[fold] = (
            totals
        )

        per_fold_summary[str(fold)] = {
            "cells":
                int(
                    observed.n_obs
                ),

            "observed_genes":
                int(
                    observed.n_vars
                ),

            "heldout_genes":
                int(
                    heldout.n_vars
                ),

            "blank_observed_cells":
                int(
                    np.sum(
                        totals == 0
                    )
                ),

            "nonfinite_library_cells":
                int(
                    np.sum(
                        ~finite
                    )
                ),
        }

        print()
        print(
            f"{donor} | Fold {fold}"
        )

        print(
            "Observed shape:",
            observed.shape,
        )

        print(
            "Held-out shape:",
            heldout.shape,
        )

        print(
            "Blank observed cells:",
            int(
                np.sum(
                    totals == 0
                )
            ),
        )

        print(
            "Library-size percentiles:",
            np.percentile(
                totals,
                [
                    0,
                    1,
                    25,
                    50,
                    75,
                    99,
                    100,
                ],
            ).tolist(),
        )

    if (
        canonical_cells is None
        or
        common_keep is None
    ):

        raise RuntimeError(
            "No fold data loaded."
        )

    # ============================================================
    # Common mask
    # ============================================================

    retained_cells = (
        canonical_cells[
            common_keep
        ]
    )

    removed_cells = (
        canonical_cells[
            ~common_keep
        ]
    )

    audit = pd.DataFrame({
        "cell_id":
            canonical_cells,

        "fold_1_observed_sum":
            observed_totals[1],

        "fold_2_observed_sum":
            observed_totals[2],

        "fold_3_observed_sum":
            observed_totals[3],
    })

    for fold in (1, 2, 3):

        audit[
            f"blank_in_fold_{fold}"
        ] = (
            observed_totals[fold]
            == 0
        )

    audit["common_keep"] = (
        common_keep
    )

    audit["removal_reason"] = (
        np.where(
            common_keep,
            "",
            (
                "zero observed-gene library "
                "in at least one fold"
            ),
        )
    )

    audit.to_csv(
        audit_path,
        index=False,
    )

    pd.DataFrame({
        "cell_id":
            retained_cells
    }).to_csv(
        common_mask_path,
        index=False,
    )

    audit.loc[
        ~common_keep
    ].to_csv(
        removed_path,
        index=False,
    )

    print("\n" + "=" * 100)
    print("COMMON MASK")
    print("=" * 100)

    print(
        "Original cells:",
        len(
            canonical_cells
        ),
    )

    print(
        "Retained cells:",
        len(
            retained_cells
        ),
    )

    print(
        "Removed cells:",
        len(
            removed_cells
        ),
    )

    # ============================================================
    # Rewrite all six active H5AD files
    # ============================================================

    filter_metadata = {
        "donor":
            donor,

        "filter_name":
            "common_nonempty_observed_cells",

        "definition":
            (
                "Cell has finite observed-gene "
                "library size > 0 in Fold 1, "
                "Fold 2, and Fold 3."
            ),

        "heldout_genes_used_for_filtering":
            False,

        "original_cells":
            int(
                len(
                    canonical_cells
                )
            ),

        "retained_cells":
            int(
                len(
                    retained_cells
                )
            ),

        "removed_cells":
            int(
                len(
                    removed_cells
                )
            ),

        "created_at":
            datetime.now().isoformat(),
    }

    print("\n" + "=" * 100)
    print("REWRITE FILTERED H5AD FILES")
    print("=" * 100)

    for fold in (1, 2, 3):

        observed_filtered = (
            observed_data[fold][
                common_keep
            ].copy()
        )

        heldout_filtered = (
            heldout_data[fold][
                common_keep
            ].copy()
        )

        # --------------------------------------------------------
        # Validate common cell universe
        # --------------------------------------------------------

        if not np.array_equal(
            observed_filtered
            .obs_names
            .to_numpy(),

            retained_cells,
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "retained observed cell "
                "order incorrect."
            )

        if not np.array_equal(
            observed_filtered
            .obs_names
            .to_numpy(),

            heldout_filtered
            .obs_names
            .to_numpy(),
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "filtered observed/truth "
                "order differs."
            )

        filtered_totals = row_sums(
            observed_filtered.X
        )

        if np.any(
            ~np.isfinite(
                filtered_totals
            )
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "nonfinite libraries remain."
            )

        if np.any(
            filtered_totals <= 0
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "blank cells remain."
            )

        observed_filtered.uns[
            "common_nonempty_filter"
        ] = filter_metadata

        heldout_filtered.uns[
            "common_nonempty_filter"
        ] = filter_metadata

        observed_destination = (
            active_path(
                fold,
                "observed",
            )
        )

        heldout_destination = (
            active_path(
                fold,
                "heldout",
            )
        )

        print(
            f"Fold {fold}: "
            f"{observed_filtered.shape} / "
            f"{heldout_filtered.shape}"
        )

        write_atomic(
            observed_filtered,
            observed_destination,
        )

        write_atomic(
            heldout_filtered,
            heldout_destination,
        )

    # ============================================================
    # Final disk validation
    # ============================================================

    print("\n" + "=" * 100)
    print("FINAL ON-DISK VALIDATION")
    print("=" * 100)

    final_summary = {}

    for fold in (1, 2, 3):

        observed = ad.read_h5ad(
            active_path(
                fold,
                "observed",
            )
        )

        heldout = ad.read_h5ad(
            active_path(
                fold,
                "heldout",
            )
        )

        totals = row_sums(
            observed.X
        )

        blank_final = int(
            np.sum(
                totals <= 0
            )
        )

        if (
            observed.n_obs
            != len(
                retained_cells
            )
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "final observed cell "
                "count incorrect."
            )

        if (
            heldout.n_obs
            != len(
                retained_cells
            )
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "final heldout cell "
                "count incorrect."
            )

        if blank_final != 0:

            raise ValueError(
                f"{donor} Fold {fold}: "
                f"{blank_final} blank cells remain."
            )

        if not np.array_equal(
            observed
            .obs_names
            .to_numpy(),

            heldout
            .obs_names
            .to_numpy(),
        ):

            raise ValueError(
                f"{donor} Fold {fold}: "
                "final observed/truth "
                "cell order differs."
            )

        final_summary[str(fold)] = {
            "cells":
                int(
                    observed.n_obs
                ),

            "observed_genes":
                int(
                    observed.n_vars
                ),

            "heldout_genes":
                int(
                    heldout.n_vars
                ),

            "remaining_blank_observed_cells":
                blank_final,
        }

        print(
            f"Fold {fold}: "
            f"observed={observed.shape}, "
            f"heldout={heldout.shape}, "
            f"blank={blank_final}"
        )

    # ============================================================
    # Manifest
    # ============================================================

    manifest = {
        **filter_metadata,

        "fold_directory":
            str(
                fold_dir
            ),

        "backup_directory":
            str(
                backup_dir
            ),

        "per_fold_before_filtering":
            per_fold_summary,

        "per_fold_after_filtering":
            final_summary,

        "common_mask":
            str(
                common_mask_path
            ),

        "removed_cells_file":
            str(
                removed_path
            ),

        "audit_file":
            str(
                audit_path
            ),
    }

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n"
    )

    print("\n" + "=" * 100)
    print("SUCCESS")
    print("=" * 100)

    print(
        "Donor:",
        donor,
    )

    print(
        "Common retained cells:",
        len(
            retained_cells
        ),
    )

    print(
        "Removed cells:",
        len(
            removed_cells
        ),
    )

    print(
        "Manifest:",
        manifest_path,
    )

    print(
        "All three observed folds now "
        "contain zero blank cells."
    )

    print(
        "\nFINAL STATUS: PASS"
    )


if __name__ == "__main__":
    main()
