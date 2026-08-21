#!/usr/bin/env python3

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

RAW_FILE = Path(
    "data/processed/imputation_beta/Br8667/"
    "spatial_data_xenium_Br8667_vista.h5ad"
)

QC_FILE = Path(
    "data/processed/imputation_beta/Br8667/"
    "spatial_data_xenium_Br8667_vista_qc.h5ad"
)

OUT_DIR = Path(
    "outputs/data_checks/imputation_beta/Br8667/"
    "xenium_qc_comparison"
)

EXPECTED_GENES = 300
ATOL = 1e-8
RTOL = 1e-5


def print_section(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def save_list(values, path, column):
    pd.DataFrame({column: list(values)}).to_csv(path, index=False)


def compare_matrix(a, b):
    if a.shape != b.shape:
        return {
            "same_shape": False,
            "identical": False,
            "max_abs_difference": np.nan,
            "mean_abs_difference": np.nan,
            "different_values": np.nan,
        }

    if sparse.issparse(a) or sparse.issparse(b):
        a = a.tocsr() if sparse.issparse(a) else sparse.csr_matrix(a)
        b = b.tocsr() if sparse.issparse(b) else sparse.csr_matrix(b)

        diff = (a.astype(float) - b.astype(float)).tocsr()
        max_abs = float(np.abs(diff.data).max()) if diff.nnz else 0.0

        union = ((a != 0) + (b != 0)).astype(bool).tocsr()
        rows, cols = union.nonzero()

        if len(rows) == 0:
            return {
                "same_shape": True,
                "identical": True,
                "max_abs_difference": 0.0,
                "mean_abs_difference": 0.0,
                "different_values": 0,
            }

        a_values = np.asarray(a[rows, cols]).ravel()
        b_values = np.asarray(b[rows, cols]).ravel()

        close = np.isclose(
            a_values,
            b_values,
            atol=ATOL,
            rtol=RTOL,
            equal_nan=True,
        )

        return {
            "same_shape": True,
            "identical": bool(close.all()),
            "max_abs_difference": max_abs,
            "mean_abs_difference": float(
                np.abs(a_values - b_values).mean()
            ),
            "different_values": int((~close).sum()),
        }

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    close = np.isclose(
        a,
        b,
        atol=ATOL,
        rtol=RTOL,
        equal_nan=True,
    )

    difference = np.abs(a - b)

    return {
        "same_shape": True,
        "identical": bool(close.all()),
        "max_abs_difference": float(np.nanmax(difference)),
        "mean_abs_difference": float(np.nanmean(difference)),
        "different_values": int((~close).sum()),
    }


def find_annotation_columns(columns):
    cell_type_columns = []
    layer_columns = []

    for column in columns:
        lower = column.lower()

        if any(
            term in lower
            for term in [
                "cell_type",
                "celltype",
                "cell type",
                "cell_annotation",
                "annotation",
                "broad_type",
            ]
        ):
            cell_type_columns.append(column)

        if any(
            term in lower
            for term in [
                "layer",
                "layerwise",
                "cortical_layer",
                "spatial_domain",
            ]
        ):
            layer_columns.append(column)

    return sorted(set(cell_type_columns)), sorted(set(layer_columns))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print_section("INPUT FILES")

    print("Raw Xenium file:")
    print(RAW_FILE.resolve())

    print("\nQC Xenium file:")
    print(QC_FILE.resolve())

    if not RAW_FILE.exists():
        raise FileNotFoundError(f"Raw file not found: {RAW_FILE}")

    if not QC_FILE.exists():
        raise FileNotFoundError(f"QC file not found: {QC_FILE}")

    print_section("LOAD ANNDATA")

    print("Loading raw file...")
    raw = ad.read_h5ad(RAW_FILE)

    print("Loading QC file...")
    qc = ad.read_h5ad(QC_FILE)

    raw_cells = pd.Index(raw.obs_names.astype(str))
    qc_cells = pd.Index(qc.obs_names.astype(str))

    raw_genes = pd.Index(raw.var_names.astype(str))
    qc_genes = pd.Index(qc.var_names.astype(str))

    common_cells = raw_cells.intersection(qc_cells, sort=False)
    removed_cells = raw_cells.difference(qc_cells, sort=False)
    extra_cells = qc_cells.difference(raw_cells, sort=False)

    common_genes = raw_genes.intersection(qc_genes, sort=False)
    missing_genes = raw_genes.difference(qc_genes, sort=False)
    extra_genes = qc_genes.difference(raw_genes, sort=False)

    same_gene_set = len(missing_genes) == 0 and len(extra_genes) == 0
    same_gene_order = raw_genes.equals(qc_genes)
    qc_cells_subset = len(extra_cells) == 0

    print_section("1. DIMENSIONS")

    print(f"Raw shape: {raw.n_obs:,} cells x {raw.n_vars:,} genes")
    print(f"QC shape:  {qc.n_obs:,} cells x {qc.n_vars:,} genes")

    print(f"\nCells removed by QC: {len(removed_cells):,}")
    print(f"Unexpected extra QC cells: {len(extra_cells):,}")

    print_section("2. GENE CHECK")

    print(f"Expected genes: {EXPECTED_GENES}")
    print(f"Raw genes:      {raw.n_vars}")
    print(f"QC genes:       {qc.n_vars}")

    print(f"\nRaw has {EXPECTED_GENES} genes: {raw.n_vars == EXPECTED_GENES}")
    print(f"QC has {EXPECTED_GENES} genes:  {qc.n_vars == EXPECTED_GENES}")
    print(f"Same gene set:               {same_gene_set}")
    print(f"Same gene order:             {same_gene_order}")

    print(f"\nCommon genes:     {len(common_genes)}")
    print(f"Missing from QC:  {len(missing_genes)}")
    print(f"Extra in QC:      {len(extra_genes)}")

    if len(missing_genes):
        print("\nGenes missing from QC:")
        print(list(missing_genes[:30]))

    if len(extra_genes):
        print("\nGenes extra in QC:")
        print(list(extra_genes[:30]))

    save_list(raw_genes, OUT_DIR / "raw_genes.csv", "gene")
    save_list(qc_genes, OUT_DIR / "qc_genes.csv", "gene")
    save_list(
        missing_genes,
        OUT_DIR / "genes_missing_from_qc.csv",
        "gene",
    )
    save_list(
        extra_genes,
        OUT_DIR / "genes_extra_in_qc.csv",
        "gene",
    )

    print_section("3. CELL CHECK")

    print(f"Raw cells:              {len(raw_cells):,}")
    print(f"QC cells:               {len(qc_cells):,}")
    print(f"Common cells:           {len(common_cells):,}")
    print(f"Cells removed by QC:    {len(removed_cells):,}")
    print(f"Unexpected QC cells:    {len(extra_cells):,}")
    print(f"QC subset of raw cells: {qc_cells_subset}")

    save_list(
        removed_cells,
        OUT_DIR / "cells_removed_by_qc.csv",
        "cell_id",
    )
    save_list(
        extra_cells,
        OUT_DIR / "cells_extra_in_qc.csv",
        "cell_id",
    )

    print_section("4. EXPRESSION CHECK")

    expression_result = {
        "identical": False,
        "max_abs_difference": np.nan,
        "mean_abs_difference": np.nan,
        "different_values": np.nan,
    }

    if len(common_cells) and len(common_genes):
        gene_order = raw_genes[raw_genes.isin(common_genes)]

        raw_common = raw[common_cells, gene_order].X
        qc_common = qc[common_cells, gene_order].X

        expression_result = compare_matrix(raw_common, qc_common)

        print(f"Shared cells compared: {len(common_cells):,}")
        print(f"Shared genes compared: {len(gene_order):,}")
        print(
            "Expression values preserved: "
            f"{expression_result['identical']}"
        )
        print(
            "Different values: "
            f"{expression_result['different_values']}"
        )
        print(
            "Maximum absolute difference: "
            f"{expression_result['max_abs_difference']}"
        )
        print(
            "Mean absolute difference: "
            f"{expression_result['mean_abs_difference']}"
        )
    else:
        print("No common cells or genes available for comparison.")

    print_section("5. SPATIAL COORDINATES")

    coordinates_available = False
    coordinates_preserved = False
    coordinate_source = None

    if "spatial" in raw.obsm and "spatial" in qc.obsm:
        coordinate_source = 'obsm["spatial"]'
        coordinates_available = True

        raw_coordinates = np.asarray(
            raw[common_cells].obsm["spatial"]
        )
        qc_coordinates = np.asarray(
            qc[common_cells].obsm["spatial"]
        )

        coordinates_preserved = bool(
            np.allclose(
                raw_coordinates,
                qc_coordinates,
                atol=ATOL,
                rtol=0,
                equal_nan=True,
            )
        )

    else:
        coordinate_pairs = [
            ("x_centroid", "y_centroid"),
            ("x_location", "y_location"),
            ("x", "y"),
            ("spatial_x", "spatial_y"),
        ]

        for x_column, y_column in coordinate_pairs:
            if (
                x_column in raw.obs.columns
                and y_column in raw.obs.columns
                and x_column in qc.obs.columns
                and y_column in qc.obs.columns
            ):
                coordinate_source = (
                    f"obs[['{x_column}', '{y_column}']]"
                )
                coordinates_available = True

                raw_coordinates = raw.obs.loc[
                    common_cells,
                    [x_column, y_column],
                ].to_numpy(dtype=float)

                qc_coordinates = qc.obs.loc[
                    common_cells,
                    [x_column, y_column],
                ].to_numpy(dtype=float)

                coordinates_preserved = bool(
                    np.allclose(
                        raw_coordinates,
                        qc_coordinates,
                        atol=ATOL,
                        rtol=0,
                        equal_nan=True,
                    )
                )
                break

    print(f"Coordinate source:             {coordinate_source}")
    print(f"Coordinates available:        {coordinates_available}")
    print(f"Coordinates preserved:        {coordinates_preserved}")

    print_section("6. OBS METADATA")

    raw_obs_columns = pd.Index(raw.obs.columns.astype(str))
    qc_obs_columns = pd.Index(qc.obs.columns.astype(str))

    obs_added = qc_obs_columns.difference(
        raw_obs_columns,
        sort=False,
    )

    obs_removed = raw_obs_columns.difference(
        qc_obs_columns,
        sort=False,
    )

    print("Raw obs columns:")
    print(list(raw_obs_columns))

    print("\nQC obs columns:")
    print(list(qc_obs_columns))

    print("\nColumns added in QC:")
    print(list(obs_added))

    print("\nColumns missing from QC:")
    print(list(obs_removed))

    save_list(
        obs_added,
        OUT_DIR / "obs_columns_added_in_qc.csv",
        "column",
    )

    save_list(
        obs_removed,
        OUT_DIR / "obs_columns_missing_from_qc.csv",
        "column",
    )

    print_section("7. CELL AND LAYER ANNOTATION")

    cell_type_columns, layer_columns = find_annotation_columns(
        qc_obs_columns
    )

    print("Candidate cell-type columns:")
    print(cell_type_columns)

    print("\nCandidate layer columns:")
    print(layer_columns)

    annotation_tables = []

    for column in cell_type_columns + layer_columns:
        counts = (
            qc.obs[column]
            .astype("string")
            .fillna("<NA>")
            .value_counts(dropna=False)
            .rename_axis("value")
            .reset_index(name="count")
        )

        counts.insert(0, "column", column)

        annotation_type = (
            "cell_type"
            if column in cell_type_columns
            else "layer"
        )

        counts.insert(0, "annotation_type", annotation_type)
        annotation_tables.append(counts)

        print(f"\n{column}:")
        print(counts.head(30).to_string(index=False))

    if annotation_tables:
        pd.concat(
            annotation_tables,
            ignore_index=True,
        ).to_csv(
            OUT_DIR / "annotation_value_counts.csv",
            index=False,
        )

    print_section("8. ANNDATA STRUCTURE")

    structures = {
        "layers": (raw.layers.keys(), qc.layers.keys()),
        "obsm": (raw.obsm.keys(), qc.obsm.keys()),
        "obsp": (raw.obsp.keys(), qc.obsp.keys()),
        "varm": (raw.varm.keys(), qc.varm.keys()),
        "varp": (raw.varp.keys(), qc.varp.keys()),
        "uns": (raw.uns.keys(), qc.uns.keys()),
    }

    structure_rows = []

    for structure_name, (raw_keys, qc_keys) in structures.items():
        raw_keys = set(map(str, raw_keys))
        qc_keys = set(map(str, qc_keys))

        added = sorted(qc_keys - raw_keys)
        removed = sorted(raw_keys - qc_keys)

        print(f"\n{structure_name}")
        print(f"Raw:             {sorted(raw_keys)}")
        print(f"QC:              {sorted(qc_keys)}")
        print(f"Added in QC:     {added}")
        print(f"Missing from QC: {removed}")

        for key in sorted(raw_keys):
            structure_rows.append(
                {
                    "structure": structure_name,
                    "status": "raw",
                    "key": key,
                }
            )

        for key in sorted(qc_keys):
            structure_rows.append(
                {
                    "structure": structure_name,
                    "status": "qc",
                    "key": key,
                }
            )

        for key in added:
            structure_rows.append(
                {
                    "structure": structure_name,
                    "status": "added_in_qc",
                    "key": key,
                }
            )

        for key in removed:
            structure_rows.append(
                {
                    "structure": structure_name,
                    "status": "missing_from_qc",
                    "key": key,
                }
            )

    pd.DataFrame(structure_rows).to_csv(
        OUT_DIR / "anndata_structure_keys.csv",
        index=False,
    )

    print_section("9. FINAL VERDICT")

    checks = {
        "raw_has_300_genes": raw.n_vars == EXPECTED_GENES,
        "qc_has_300_genes": qc.n_vars == EXPECTED_GENES,
        "same_gene_set": same_gene_set,
        "same_gene_order": same_gene_order,
        "qc_cells_subset_of_raw": qc_cells_subset,
        "expression_preserved": expression_result["identical"],
        "cell_type_annotation_found": len(cell_type_columns) > 0,
        "layer_annotation_found": len(layer_columns) > 0,
    }

    if coordinates_available:
        checks["coordinates_preserved"] = coordinates_preserved

    for name, passed in checks.items():
        status = "PASS" if passed else "FAIL"
        print(f"{status}: {name}")

    critical_checks = [
        checks["qc_has_300_genes"],
        checks["same_gene_set"],
        checks["same_gene_order"],
        checks["qc_cells_subset_of_raw"],
        checks["expression_preserved"],
    ]

    annotation_checks = [
        checks["cell_type_annotation_found"],
        checks["layer_annotation_found"],
    ]

    critical_pass = all(critical_checks)
    annotation_pass = all(annotation_checks)
    overall_pass = critical_pass and annotation_pass

    print("\nCritical data integrity:")
    print("PASS" if critical_pass else "FAIL")

    print("\nAnnotation check:")
    print("PASS" if annotation_pass else "FAIL")

    print("\nOVERALL:")
    if overall_pass:
        print(
            "PASS — QC Xenium retains all 300 genes, preserves "
            "expression for shared cells, and includes cell and "
            "layer annotation."
        )
    else:
        print(
            "FAIL — Review the detailed results before starting "
            "imputation."
        )

    summary = pd.DataFrame(
        [
            ["raw_cells", raw.n_obs],
            ["qc_cells", qc.n_obs],
            ["cells_removed", len(removed_cells)],
            ["raw_genes", raw.n_vars],
            ["qc_genes", qc.n_vars],
            ["missing_genes", len(missing_genes)],
            ["extra_genes", len(extra_genes)],
            ["same_gene_set", same_gene_set],
            ["same_gene_order", same_gene_order],
            ["expression_preserved", expression_result["identical"]],
            ["coordinates_available", coordinates_available],
            ["coordinates_preserved", coordinates_preserved],
            ["cell_type_columns", ";".join(cell_type_columns)],
            ["layer_columns", ";".join(layer_columns)],
            ["critical_pass", critical_pass],
            ["annotation_pass", annotation_pass],
            ["overall_pass", overall_pass],
        ],
        columns=["metric", "value"],
    )

    summary.to_csv(
        OUT_DIR / "comparison_summary.csv",
        index=False,
    )

    print_section("OUTPUT DIRECTORY")
    print(OUT_DIR.resolve())

    if not overall_pass:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
