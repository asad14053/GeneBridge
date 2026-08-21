#!/usr/bin/env python

from pathlib import Path
import argparse

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)
#########################just full snRNA NTC######################
INPUT_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

########################If consider SCZ overllaped genes############
#INPUT_H5AD = (
#    PROJECT_ROOT
#    / "data"
#    / "processed"
#    / "snrnaseq"
#    / "sce_DLPFC_annotated"
#    / "huuki_snrna_reference_baituk_overlap.h5ad"
#)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--brain-id",
        "--brain_id",
        dest="brain_id",
        default="Br8667",
    )

    parser.add_argument(
        "--top-genes",
        "--top_genes",
        dest="top_genes",
        default="full",
    )

    parser.add_argument(
        "--compression",
        default="lzf",
        choices=["lzf", "gzip", "none"],
    )

    return parser.parse_args()


def main():

    args = parse_args()

    # -------------------------------------------------------------------------
    # Brain ID
    # -------------------------------------------------------------------------

    brain_clean = str(args.brain_id).strip()

    if brain_clean.lower().startswith("br"):
        brain_clean = brain_clean[2:]

    brain_label = f"Br{brain_clean}"


    # -------------------------------------------------------------------------
    # Gene mode
    # -------------------------------------------------------------------------

    top_genes_value = str(args.top_genes).lower()

    if top_genes_value in {
        "full",
        "all",
        "none",
        "inf",
        "infinite",
    }:

        top_n = None
        gene_mode = "full"

    else:

        top_n = int(top_genes_value)

        if top_n <= 0:
            raise ValueError(
                "--top-genes must be a positive integer or 'full'"
            )

        gene_mode = f"top{top_n}"


    # -------------------------------------------------------------------------
    # Output paths
    # -------------------------------------------------------------------------

    OUT_DIR = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "imputation_beta"
        / brain_label
    )

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUT_H5AD = (
        OUT_DIR
        / f"seq_data_huuki_snrna_{brain_label}_{gene_mode}.h5ad"
    )

    OUT_GENES = (
        OUT_DIR
        / f"seq_data_huuki_snrna_{brain_label}_{gene_mode}_genes.csv"
    )


    # -------------------------------------------------------------------------
    # Print configuration
    # -------------------------------------------------------------------------

    print("=" * 80)
    print("Huuki donor-specific snRNA-seq H5AD creation")
    print("=" * 80)

    print("Input:", INPUT_H5AD)
    print("Brain ID:", brain_label)
    print("Gene mode:", gene_mode)
    print("Output:", OUT_H5AD)
    print("Gene list:", OUT_GENES)


    # -------------------------------------------------------------------------
    # Check input
    # -------------------------------------------------------------------------

    if not INPUT_H5AD.exists():
        raise FileNotFoundError(
            f"Input H5AD not found:\n{INPUT_H5AD}"
        )


    # -------------------------------------------------------------------------
    # Read full Huuki H5AD in backed mode
    # -------------------------------------------------------------------------

    print("\nReading full Huuki reference...")

    adata = ad.read_h5ad(
        INPUT_H5AD,
        backed="r",
    )

    print(adata)

    print(
        "Full shape [cells x genes]:",
        adata.shape,
    )

    print(
        "Layers:",
        list(adata.layers.keys()),
    )


    # -------------------------------------------------------------------------
    # Check BrNum
    # -------------------------------------------------------------------------

    if "BrNum" not in adata.obs.columns:

        print(
            "\nAvailable obs columns:"
        )

        print(
            list(adata.obs.columns)
        )

        adata.file.close()

        raise KeyError(
            "BrNum not found in adata.obs"
        )


    # -------------------------------------------------------------------------
    # Select Br8667 cells
    # -------------------------------------------------------------------------

    br_values = (
        adata.obs["BrNum"]
        .astype(str)
        .str.replace(
            r"^Br",
            "",
            regex=True,
        )
    )

    keep_cells = (
        br_values
        == brain_clean
    )

    cell_indices = np.flatnonzero(
        keep_cells.to_numpy()
    )

    print(
        f"\nCells found for {brain_label}:",
        len(cell_indices),
    )


    if len(cell_indices) == 0:

        print(
            "Available BrNum values:"
        )

        print(
            sorted(
                adata.obs["BrNum"]
                .astype(str)
                .unique()
                .tolist()
            )
        )

        adata.file.close()

        raise ValueError(
            f"No cells found for {brain_label}"
        )


    # -------------------------------------------------------------------------
    # Select genes
    # -------------------------------------------------------------------------

    if top_n is None:

        print(
            "\nKeeping FULL gene set."
        )

        gene_indices = np.arange(
            adata.n_vars
        )

    else:

        if (
            "binomial_deviance"
            not in adata.var.columns
        ):

            adata.file.close()

            raise KeyError(
                "binomial_deviance not found in adata.var"
            )

        print(
            "\nSelecting top genes by binomial_deviance:"
        )

        print(
            top_n
        )

        deviance = pd.to_numeric(
            adata.var["binomial_deviance"],
            errors="coerce",
        ).to_numpy()

        valid = np.flatnonzero(
            np.isfinite(deviance)
        )

        ranking = valid[
            np.argsort(
                -deviance[valid],
                kind="stable",
            )
        ]

        gene_indices = ranking[
            :min(
                top_n,
                len(ranking),
            )
        ]


    print(
        "Genes selected:",
        len(gene_indices),
    )


    # -------------------------------------------------------------------------
    # Load ONLY selected donor and genes into memory
    # -------------------------------------------------------------------------

    print(
        "\nLoading donor cells into memory..."
    )

    # Subset only cells while backed.
    # HDF5/h5py does not allow fancy indexing on both axes at once.
    adata_b = (
        adata[
            cell_indices,
            :
        ]
        .to_memory()
    )

    adata.file.close()

    print(
        "Donor subset loaded [cells x genes]:",
        adata_b.shape,
    )

    # Subset genes only after the donor subset is in memory.
    if top_n is not None:

        print(
            "\nSubsetting selected genes in memory..."
        )

        adata_b = (
            adata_b[
                :,
                gene_indices
            ]
            .copy()
        )

    print(
        "Final subset shape [cells x genes]:",
        adata_b.shape,
    )


    print(
        "Subset shape [cells x genes]:",
        adata_b.shape,
    )


    # -------------------------------------------------------------------------
    # Match previous R output:
    #
    # zellkonverter::writeH5AD(
    #     sce_b,
    #     X_name = "logcounts"
    # )
    #
    # Therefore final X should be logcounts.
    # Preserve counts in layers["counts"].
    # -------------------------------------------------------------------------

    print(
        "\nPreparing X = logcounts..."
    )


    if "logcounts" not in adata_b.layers:

        raise KeyError(
            "logcounts layer not found in merged Huuki H5AD"
        )


    # Preserve current X as counts if counts layer is absent
    if "counts" not in adata_b.layers:

        adata_b.layers["counts"] = (
            adata_b.X.copy()
        )


    # Put logcounts into X
    adata_b.X = (
        adata_b.layers[
            "logcounts"
        ].copy()
    )


    adata_b.obs_names_make_unique()

    adata_b.var_names_make_unique()


    print(
        "Final X: logcounts"
    )

    print(
        "Final layers:",
        list(adata_b.layers.keys()),
    )


    # -------------------------------------------------------------------------
    # Save gene list
    # -------------------------------------------------------------------------

    gene_df = pd.DataFrame(
        {
            "gene":
                adata_b.var_names
                .astype(str)
                .tolist()
        }
    )


    for column in [
        "gene_id",
        "gene_name",
        "binomial_deviance",
    ]:

        if column in adata_b.var.columns:

            gene_df[column] = (
                adata_b.var[
                    column
                ].to_numpy()
            )


    gene_df.to_csv(
        OUT_GENES,
        index=False,
    )


    # -------------------------------------------------------------------------
    # Remove existing H5AD
    # -------------------------------------------------------------------------

    if OUT_H5AD.exists():

        print(
            "\nRemoving existing output:"
        )

        print(
            OUT_H5AD
        )

        OUT_H5AD.unlink()


    # -------------------------------------------------------------------------
    # Write H5AD
    # -------------------------------------------------------------------------

    compression = (
        None
        if args.compression == "none"
        else args.compression
    )


    print(
        "\nWriting H5AD..."
    )

    adata_b.write_h5ad(
        OUT_H5AD,
        compression=compression,
    )


    # -------------------------------------------------------------------------
    # Verify saved output
    # -------------------------------------------------------------------------

    print(
        "\nVerifying saved H5AD..."
    )

    check = ad.read_h5ad(
        OUT_H5AD,
        backed="r",
    )

    print(check)

    print(
        "Saved shape:",
        check.shape,
    )

    print(
        "Saved layers:",
        list(check.layers.keys()),
    )

    check.file.close()


    # -------------------------------------------------------------------------
    # Done
    # -------------------------------------------------------------------------

    print(
        "\n============================================================"
    )

    print(
        "DONE"
    )

    print(
        "============================================================"
    )

    print(
        "Brain:",
        brain_label,
    )

    print(
        "Cells:",
        adata_b.n_obs,
    )

    print(
        "Genes:",
        adata_b.n_vars,
    )

    print(
        "Output H5AD:"
    )

    print(
        OUT_H5AD
    )

    print(
        "\nSaved gene list:"
    )

    print(
        OUT_GENES
    )


if __name__ == "__main__":
    main()
