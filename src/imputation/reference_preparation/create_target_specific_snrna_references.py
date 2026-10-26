#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

DEFAULT_REFERENCE = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated/"
      "huuki_snrna_reference_full_allgenes.h5ad"
)

DEFAULT_OUTPUT_ROOT = (
    PROJECT_ROOT
    / "data/processed/imputation_full/"
      "snrnaseq_target_references"
)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--donor",
        required=True,
        choices=["Br6432", "Br8667"],
    )

    parser.add_argument(
        "--reference",
        type=Path,
        default=DEFAULT_REFERENCE,
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
    )

    parser.add_argument(
        "--donor-col",
        default="BrNum",
    )

    return parser.parse_args()


def validate_gene_universe(
    label: str,
    adata: ad.AnnData,
    canonical_genes: pd.Index,
):
    genes = pd.Index(
        adata.var_names.astype(str)
    )

    if not genes.is_unique:
        raise ValueError(
            f"{label}: duplicate gene names detected."
        )

    if len(genes) != len(canonical_genes):
        raise ValueError(
            f"{label}: expected {len(canonical_genes)} genes, "
            f"found {len(genes)}."
        )

    if not genes.equals(canonical_genes):
        raise ValueError(
            f"{label}: gene identities/order differ "
            "from canonical full reference."
        )


def main():

    args = parse_args()

    donor = args.donor
    reference_path = args.reference

    output_dir = (
        args.output_root
        / donor
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ============================================================
    # CANONICAL FILENAMES
    # ============================================================

    matched_path = (
        output_dir
        / f"seq_data_huuki_snrna_{donor}.h5ad"
    )

    excluded_path = (
        output_dir
        / f"seq_data_huuki_snrna_9brains_excluding_{donor}.h5ad"
    )

    manifest_path = (
        output_dir
        / "reference_preparation_manifest.json"
    )

    print("=" * 100)
    print("CREATE TARGET-SPECIFIC snRNA REFERENCES")
    print("=" * 100)

    print("Target donor :", donor)
    print("Source       :", reference_path)
    print("Output dir   :", output_dir)

    if not reference_path.exists():
        raise FileNotFoundError(
            reference_path
        )

    # ============================================================
    # Load canonical 10-brain reference
    # ============================================================

    ref = ad.read_h5ad(
        reference_path,
        backed="r",
    )

    print()
    print("=" * 100)
    print("CANONICAL FULL REFERENCE")
    print("=" * 100)

    print("Shape:", ref.shape)

    if args.donor_col not in ref.obs.columns:
        raise ValueError(
            f"Missing obs column: {args.donor_col}"
        )

    if not ref.obs_names.is_unique:
        raise ValueError(
            "Reference cell IDs are not unique."
        )

    if not ref.var_names.is_unique:
        raise ValueError(
            "Reference gene names are not unique."
        )

    brain = (
        ref.obs[args.donor_col]
        .astype(str)
        .str.strip()
    )

    canonical_genes = pd.Index(
        ref.var_names.astype(str)
    )

    brain_counts = (
        brain
        .value_counts()
        .sort_index()
    )

    brains = sorted(
        brain.unique()
    )

    print("Brains:", len(brains))
    print("Genes :", len(canonical_genes))

    print()
    print("Cells per brain:")
    print(
        brain_counts.to_string()
    )

    if len(brains) != 10:
        raise ValueError(
            f"Expected 10 brains, found {len(brains)}."
        )

    if donor not in brains:
        raise ValueError(
            f"{donor} is not present in the "
            "canonical 10-brain reference."
        )

    # ============================================================
    # Build masks
    # ============================================================

    donor_array = brain.to_numpy()

    matched_mask = (
        donor_array == donor
    )

    excluded_mask = (
        donor_array != donor
    )

    n_matched = int(
        matched_mask.sum()
    )

    n_excluded = int(
        excluded_mask.sum()
    )

    print()
    print("=" * 100)
    print("SUBSET DESIGN")
    print("=" * 100)

    print(
        "Full reference cells:",
        ref.n_obs,
    )

    print(
        f"{donor}-only cells:",
        n_matched,
    )

    print(
        f"9-brain excluding {donor}:",
        n_excluded,
    )

    print(
        "Partition check:",
        n_matched + n_excluded,
    )

    if (
        n_matched
        + n_excluded
        != ref.n_obs
    ):
        raise ValueError(
            "Matched + excluded cell counts "
            "do not equal full reference."
        )

    if n_matched == 0:
        raise ValueError(
            f"No cells found for {donor}."
        )

    # ============================================================
    # Ex-5: matched donor only
    # ============================================================

    print()
    print("=" * 100)
    print("EX-5: MATCHED SINGLE-DONOR REFERENCE")
    print("=" * 100)

    matched = (
        ref[
            np.flatnonzero(
                matched_mask
            ),
            :
        ]
        .to_memory()
    )

    matched.obs_names = (
        matched.obs_names.astype(str)
    )

    matched.var_names = (
        matched.var_names.astype(str)
    )

    validate_gene_universe(
        "Matched reference",
        matched,
        canonical_genes,
    )

    matched_brains = sorted(
        matched.obs[
            args.donor_col
        ]
        .astype(str)
        .str.strip()
        .unique()
    )

    if matched_brains != [donor]:
        raise ValueError(
            f"Matched reference contains unexpected brains: "
            f"{matched_brains}"
        )

    matched.uns[
        "target_specific_reference"
    ] = {
        "target_donor":
            donor,

        "experiment":
            "ex5",

        "reference_role":
            "matched_single_donor",

        "source_reference":
            str(
                reference_path.resolve()
            ),

        "created_at":
            datetime.now().isoformat(),
    }

    print(
        "Shape:",
        matched.shape,
    )

    print(
        "Brains:",
        matched_brains,
    )

    print(
        "Writing:",
        matched_path,
    )

    matched.write_h5ad(
        matched_path,
        compression="gzip",
    )

    del matched

    # ============================================================
    # Ex-5.3: 9 brains excluding target
    # ============================================================

    print()
    print("=" * 100)
    print("EX-5.3: NINE-BRAIN TARGET-EXCLUDED REFERENCE")
    print("=" * 100)

    excluded = (
        ref[
            np.flatnonzero(
                excluded_mask
            ),
            :
        ]
        .to_memory()
    )

    excluded.obs_names = (
        excluded.obs_names.astype(str)
    )

    excluded.var_names = (
        excluded.var_names.astype(str)
    )

    validate_gene_universe(
        "Excluded reference",
        excluded,
        canonical_genes,
    )

    excluded_brains = sorted(
        excluded.obs[
            args.donor_col
        ]
        .astype(str)
        .str.strip()
        .unique()
    )

    if donor in excluded_brains:
        raise ValueError(
            f"{donor} is still present in "
            "the target-excluded reference."
        )

    if len(excluded_brains) != 9:
        raise ValueError(
            f"Expected 9 brains after excluding {donor}, "
            f"found {len(excluded_brains)}."
        )

    excluded.uns[
        "target_specific_reference"
    ] = {
        "target_donor":
            donor,

        "experiment":
            "ex5_3",

        "reference_role":
            "nine_brains_excluding_target",

        "source_reference":
            str(
                reference_path.resolve()
            ),

        "created_at":
            datetime.now().isoformat(),
    }

    print(
        "Shape:",
        excluded.shape,
    )

    print(
        "Brains:",
        excluded_brains,
    )

    print(
        "Writing:",
        excluded_path,
    )

    excluded.write_h5ad(
        excluded_path,
        compression="gzip",
    )

    del excluded

    ref.file.close()

    # ============================================================
    # Final on-disk validation
    # ============================================================

    print()
    print("=" * 100)
    print("FINAL ON-DISK VALIDATION")
    print("=" * 100)

    matched_check = ad.read_h5ad(
        matched_path,
        backed="r",
    )

    excluded_check = ad.read_h5ad(
        excluded_path,
        backed="r",
    )

    validate_gene_universe(
        "Matched on-disk",
        matched_check,
        canonical_genes,
    )

    validate_gene_universe(
        "Excluded on-disk",
        excluded_check,
        canonical_genes,
    )

    matched_brains_check = sorted(
        matched_check.obs[
            args.donor_col
        ]
        .astype(str)
        .str.strip()
        .unique()
    )

    excluded_brains_check = sorted(
        excluded_check.obs[
            args.donor_col
        ]
        .astype(str)
        .str.strip()
        .unique()
    )

    if matched_brains_check != [donor]:
        raise ValueError(
            "Final matched-reference validation failed."
        )

    if donor in excluded_brains_check:
        raise ValueError(
            "Target donor remains in final excluded reference."
        )

    if len(excluded_brains_check) != 9:
        raise ValueError(
            "Final excluded reference does not contain 9 brains."
        )

    if (
        matched_check.n_obs
        + excluded_check.n_obs
        != (
            n_matched
            + n_excluded
        )
    ):
        raise ValueError(
            "Final cell partition validation failed."
        )

    if not pd.Index(
        matched_check.var_names.astype(str)
    ).equals(
        pd.Index(
            excluded_check.var_names.astype(str)
        )
    ):
        raise ValueError(
            "Matched and excluded references "
            "do not have identical gene universes."
        )

    print()
    print("MATCHED REFERENCE")
    print(
        "File  :",
        matched_path,
    )
    print(
        "Shape :",
        matched_check.shape,
    )
    print(
        "Brains:",
        matched_brains_check,
    )

    print()
    print("EXCLUDED REFERENCE")
    print(
        "File  :",
        excluded_path,
    )
    print(
        "Shape :",
        excluded_check.shape,
    )
    print(
        "Brains:",
        excluded_brains_check,
    )

    matched_shape = [
        int(matched_check.n_obs),
        int(matched_check.n_vars),
    ]

    excluded_shape = [
        int(excluded_check.n_obs),
        int(excluded_check.n_vars),
    ]

    matched_check.file.close()
    excluded_check.file.close()

    # ============================================================
    # Manifest
    # ============================================================

    manifest = {
        "target_donor":
            donor,

        "source_reference":
            str(
                reference_path.resolve()
            ),

        "source_brains":
            brains,

        "source_cells":
            int(
                n_matched
                + n_excluded
            ),

        "source_genes":
            int(
                len(
                    canonical_genes
                )
            ),

        "ex5": {
            "reference":
                str(
                    matched_path.resolve()
                ),

            "shape":
                matched_shape,

            "brains": [
                donor
            ],
        },

        "ex5_1": {
            "reference":
                str(
                    reference_path.resolve()
                ),

            "brains":
                brains,
        },

        "ex5_3": {
            "reference":
                str(
                    excluded_path.resolve()
                ),

            "shape":
                excluded_shape,

            "brains":
                excluded_brains_check,
        },

        "validation": {
            "source_has_10_brains":
                True,

            "matched_has_1_brain":
                True,

            "excluded_has_9_brains":
                True,

            "target_absent_from_excluded":
                True,

            "identical_gene_universe":
                True,

            "exact_cell_partition":
                True,
        },

        "created_at":
            datetime.now().isoformat(),
    }

    manifest_path.write_text(
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

    print(
        "Target:",
        donor,
    )

    print()
    print("Ex-5 reference:")
    print(
        matched_path
    )

    print()
    print("Ex-5.1 reference:")
    print(
        reference_path
    )

    print()
    print("Ex-5.3 reference:")
    print(
        excluded_path
    )

    print()
    print(
        "Manifest:",
        manifest_path,
    )

    print()
    print(
        "FINAL STATUS: PASS"
    )


if __name__ == "__main__":
    main()
