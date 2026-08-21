#!/usr/bin/env python3

from __future__ import annotations

import json
import os
import shutil
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

FOLD_DIR = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667/gene_folds_200_100"
)

BACKUP_DIR = FOLD_DIR / "backup_before_common_nonempty_filter"

COMMON_MASK_PATH = FOLD_DIR / "common_nonempty_cells.csv"
REMOVED_PATH = FOLD_DIR / "removed_blank_cells.csv"
AUDIT_PATH = FOLD_DIR / "blank_cell_audit_all_cells.csv"
MANIFEST_PATH = FOLD_DIR / "blank_cell_filter_manifest.json"


def row_sums(matrix) -> np.ndarray:
    """Return one total per cell."""
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()

    return np.asarray(matrix).sum(axis=1)


def matrix_values(matrix) -> np.ndarray:
    """Return stored values for finite/nonnegative validation."""
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)

    return np.asarray(matrix).ravel()


def validate_expression(name: str, matrix) -> None:
    """Fail immediately for invalid model inputs."""
    values = matrix_values(matrix)

    if values.size == 0:
        return

    nonfinite = int(np.sum(~np.isfinite(values)))
    negative = int(np.sum(values < 0))

    if nonfinite:
        raise ValueError(
            f"{name}: found {nonfinite} nonfinite expression values."
        )

    if negative:
        raise ValueError(
            f"{name}: found {negative} negative expression values."
        )


def active_path(fold: int, kind: str) -> Path:
    return FOLD_DIR / f"fold_{fold}_{kind}_genes.h5ad"


def backup_path(fold: int, kind: str) -> Path:
    return BACKUP_DIR / f"fold_{fold}_{kind}_genes.h5ad"


def create_backups() -> None:
    """
    Preserve the original unfiltered fold files.

    On later executions, the mask is reconstructed from these backups rather
    than from already-filtered active files.
    """
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)

    for fold in (1, 2, 3):
        for kind in ("observed", "heldout"):
            source = active_path(fold, kind)
            backup = backup_path(fold, kind)

            if not source.exists():
                raise FileNotFoundError(f"Missing fold file: {source}")

            if backup.exists():
                print(f"Backup already exists: {backup}")
            else:
                print(f"Creating backup: {backup}")
                shutil.copy2(source, backup)


def write_atomic(adata: ad.AnnData, destination: Path) -> None:
    """Write to a temporary file and replace only after successful writing."""
    temporary = destination.with_name(
        destination.stem + ".tmp.h5ad"
    )

    if temporary.exists():
        temporary.unlink()

    adata.write_h5ad(temporary, compression="gzip")
    os.replace(temporary, destination)


def main() -> None:
    print("=" * 100)
    print("COMMON NONEMPTY-CELL FILTER FOR ALL THREE FOLDS")
    print("=" * 100)
    print("Fold directory:", FOLD_DIR)

    create_backups()

    observed_data: dict[int, ad.AnnData] = {}
    heldout_data: dict[int, ad.AnnData] = {}
    observed_totals: dict[int, np.ndarray] = {}

    canonical_cells: np.ndarray | None = None
    common_keep: np.ndarray | None = None
    per_fold_summary: dict[str, dict[str, int]] = {}

    # ------------------------------------------------------------------
    # Read the original backed-up fold files.
    # ------------------------------------------------------------------
    for fold in (1, 2, 3):
        observed = ad.read_h5ad(backup_path(fold, "observed"))
        heldout = ad.read_h5ad(backup_path(fold, "heldout"))

        observed.obs_names = observed.obs_names.astype(str)
        heldout.obs_names = heldout.obs_names.astype(str)

        if observed.n_vars != 200:
            raise ValueError(
                f"Fold {fold}: expected 200 observed genes, "
                f"found {observed.n_vars}."
            )

        if heldout.n_vars != 100:
            raise ValueError(
                f"Fold {fold}: expected 100 held-out genes, "
                f"found {heldout.n_vars}."
            )

        if not np.array_equal(
            observed.obs_names.to_numpy(),
            heldout.obs_names.to_numpy(),
        ):
            raise ValueError(
                f"Fold {fold}: observed and held-out cell order differs."
            )

        validate_expression(
            f"Fold {fold} observed",
            observed.X,
        )
        validate_expression(
            f"Fold {fold} held-out",
            heldout.X,
        )

        cells = observed.obs_names.to_numpy()
        totals = row_sums(observed.X)

        if canonical_cells is None:
            canonical_cells = cells.copy()
            common_keep = np.ones(
                observed.n_obs,
                dtype=bool,
            )
        elif not np.array_equal(canonical_cells, cells):
            raise ValueError(
                f"Fold {fold}: cell IDs/order differ from Fold 1."
            )

        finite = np.isfinite(totals)
        nonempty = totals > 0
        fold_keep = finite & nonempty

        common_keep &= fold_keep

        observed_data[fold] = observed
        heldout_data[fold] = heldout
        observed_totals[fold] = totals

        per_fold_summary[str(fold)] = {
            "cells": int(observed.n_obs),
            "observed_genes": int(observed.n_vars),
            "heldout_genes": int(heldout.n_vars),
            "blank_observed_cells": int(np.sum(totals == 0)),
            "nonfinite_library_cells": int(np.sum(~finite)),
        }

        print()
        print(f"Fold {fold}")
        print("-" * 50)
        print("Observed shape:", observed.shape)
        print("Held-out shape:", heldout.shape)
        print("Blank observed cells:", int(np.sum(totals == 0)))
        print(
            "Library-size percentiles:",
            np.percentile(
                totals,
                [0, 1, 25, 50, 75, 99, 100],
            ).tolist(),
        )

    if canonical_cells is None or common_keep is None:
        raise RuntimeError("No fold data were loaded.")

    retained_cells = canonical_cells[common_keep]
    removed_cells = canonical_cells[~common_keep]

    # ------------------------------------------------------------------
    # Audit table: filtering depends only on observed genes.
    # ------------------------------------------------------------------
    audit = pd.DataFrame(
        {
            "cell_id": canonical_cells,
            "fold_1_observed_sum": observed_totals[1],
            "fold_2_observed_sum": observed_totals[2],
            "fold_3_observed_sum": observed_totals[3],
        }
    )

    for fold in (1, 2, 3):
        audit[f"blank_in_fold_{fold}"] = (
            observed_totals[fold] == 0
        )

    audit["common_keep"] = common_keep

    audit["removal_reason"] = np.where(
        common_keep,
        "",
        "zero observed-gene library in at least one fold",
    )

    audit.to_csv(AUDIT_PATH, index=False)

    pd.DataFrame(
        {"cell_id": retained_cells}
    ).to_csv(COMMON_MASK_PATH, index=False)

    audit.loc[~common_keep].to_csv(
        REMOVED_PATH,
        index=False,
    )

    print()
    print("=" * 100)
    print("COMMON MASK")
    print("=" * 100)
    print("Original cells:", len(canonical_cells))
    print("Retained cells:", len(retained_cells))
    print("Removed cells:", len(removed_cells))
    print("Mask:", COMMON_MASK_PATH)
    print("Removed-cell audit:", REMOVED_PATH)

    # ------------------------------------------------------------------
    # Filter all observed and truth files using exactly the same mask.
    # ------------------------------------------------------------------
    filter_metadata = {
        "filter_name": "common_nonempty_observed_cells",
        "definition": (
            "Cell has finite observed-gene library size > 0 "
            "in Fold 1, Fold 2, and Fold 3."
        ),
        "heldout_genes_used_for_filtering": False,
        "original_cells": int(len(canonical_cells)),
        "retained_cells": int(len(retained_cells)),
        "removed_cells": int(len(removed_cells)),
        "created_at": datetime.now().isoformat(),
    }

    for fold in (1, 2, 3):
        observed_filtered = observed_data[fold][common_keep].copy()
        heldout_filtered = heldout_data[fold][common_keep].copy()

        if not np.array_equal(
            observed_filtered.obs_names.to_numpy(),
            retained_cells,
        ):
            raise ValueError(
                f"Fold {fold}: retained observed cell order is incorrect."
            )

        if not np.array_equal(
            observed_filtered.obs_names.to_numpy(),
            heldout_filtered.obs_names.to_numpy(),
        ):
            raise ValueError(
                f"Fold {fold}: filtered observed/truth order differs."
            )

        filtered_totals = row_sums(observed_filtered.X)

        if np.any(~np.isfinite(filtered_totals)):
            raise ValueError(
                f"Fold {fold}: nonfinite libraries remain after filtering."
            )

        if np.any(filtered_totals <= 0):
            raise ValueError(
                f"Fold {fold}: blank cells remain after filtering."
            )

        observed_filtered.uns[
            "common_nonempty_filter"
        ] = filter_metadata

        heldout_filtered.uns[
            "common_nonempty_filter"
        ] = filter_metadata

        observed_destination = active_path(fold, "observed")
        heldout_destination = active_path(fold, "heldout")

        print()
        print(f"Writing Fold {fold}:")
        print(" ", observed_destination)
        write_atomic(
            observed_filtered,
            observed_destination,
        )

        print(" ", heldout_destination)
        write_atomic(
            heldout_filtered,
            heldout_destination,
        )

    # ------------------------------------------------------------------
    # Reopen the active files and perform final on-disk validation.
    # ------------------------------------------------------------------
    print()
    print("=" * 100)
    print("FINAL ON-DISK VALIDATION")
    print("=" * 100)

    final_summary: dict[str, dict[str, int]] = {}

    for fold in (1, 2, 3):
        observed = ad.read_h5ad(
            active_path(fold, "observed")
        )
        heldout = ad.read_h5ad(
            active_path(fold, "heldout")
        )

        observed_totals_final = row_sums(observed.X)

        blank_final = int(
            np.sum(observed_totals_final <= 0)
        )

        if observed.n_obs != len(retained_cells):
            raise ValueError(
                f"Fold {fold}: expected {len(retained_cells)} cells, "
                f"found {observed.n_obs}."
            )

        if heldout.n_obs != len(retained_cells):
            raise ValueError(
                f"Fold {fold}: held-out cell count is incorrect."
            )

        if blank_final != 0:
            raise ValueError(
                f"Fold {fold}: {blank_final} blank cells remain."
            )

        if not np.array_equal(
            observed.obs_names.to_numpy(),
            heldout.obs_names.to_numpy(),
        ):
            raise ValueError(
                f"Fold {fold}: final observed/truth order differs."
            )

        final_summary[str(fold)] = {
            "cells": int(observed.n_obs),
            "observed_genes": int(observed.n_vars),
            "heldout_genes": int(heldout.n_vars),
            "remaining_blank_observed_cells": blank_final,
        }

        print(
            f"Fold {fold}: "
            f"observed={observed.shape}, "
            f"heldout={heldout.shape}, "
            f"blank cells={blank_final}"
        )


    manifest = {
        **filter_metadata,
        "fold_directory": str(FOLD_DIR),
        "backup_directory": str(BACKUP_DIR),
        "per_fold_before_filtering": per_fold_summary,
        "per_fold_after_filtering": final_summary,
        "common_mask": str(COMMON_MASK_PATH),
        "removed_cells_file": str(REMOVED_PATH),
        "audit_file": str(AUDIT_PATH),
    }

    MANIFEST_PATH.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n"
    )

    print()
    print("=" * 100)
    print("SUCCESS")
    print("=" * 100)
    print("Common retained cells:", len(retained_cells))
    print("Removed cells:", len(removed_cells))
    print("Manifest:", MANIFEST_PATH)
    print("All three observed folds now contain zero blank cells.")


if __name__ == "__main__":
    main()
