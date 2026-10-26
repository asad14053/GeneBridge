#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

INPUT = ROOT / "data/raw/nemo_mccarroll/counts"
OUTPUT = ROOT / "results/analysis/nemo_mccarroll/inspection"

HUUKI = (
    ROOT
    / "data/processed/snrnaseq"
    / "sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)


def decode(values):
    return [
        x.decode("utf-8", errors="replace")
        if isinstance(x, bytes) else str(x)
        for x in values
    ]


def stats(values):
    values = np.asarray(values)
    if values.size == 0:
        return {}

    return {
        "min": float(np.min(values)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def inspect_file(path, output, deep=False, chunk_cells=256):

    result = {
        "filename": path.name,
        "size_MB": round(path.stat().st_size / 1024**2, 2),
    }

    with h5py.File(path, "r") as f:

        print("\nROOT KEYS:", list(f.keys()))

        # Inspect full HDF5 hierarchy
        datasets = {}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                datasets[name] = {
                    "shape": list(obj.shape),
                    "dtype": str(obj.dtype),
                }

        f.visititems(visit)
        result["datasets"] = datasets
        result["root_attributes"] = {
            k: repr(v)[:300]
            for k, v in f.attrs.items()
        }

        print("\nDATASETS:")
        for name, info in datasets.items():
            print(name, info)

        if "matrix" not in f or "shape" not in f["matrix"]:
            result["error"] = "Missing matrix/shape"
            return result, None

        g = f["matrix"]

        required = [
            "data",
            "indices",
            "indptr",
            "barcodes",
            "features/name",
        ]

        missing = [key for key in required if key not in g]

        if missing:
            result["error"] = f"Missing datasets: {missing}"
            return result, None

        n_genes, n_cells = map(int, g["shape"][:])

        genes = decode(g["features/name"][:])
        barcodes = decode(g["barcodes"][:])

        indptr = g["indptr"][:]
        nnz = len(g["data"])

        result.update({
            "genes": n_genes,
            "nuclei": n_cells,
            "stored_entries": nnz,
            "data_dtype": str(g["data"].dtype),
            "gene_examples": genes[:10],
            "barcode_examples": barcodes[:10],
            "gene_duplicates": len(genes) - len(set(genes)),
            "barcode_duplicates": len(barcodes) - len(set(barcodes)),
            "feature_fields": list(g["features"].keys()),
            "gene_order_sha256": hashlib.sha256(
                "\n".join(genes).encode()
            ).hexdigest(),
            "matrix_shape_consistent": (
                len(genes) == n_genes
                and len(barcodes) == n_cells
            ),
            "data_indices_lengths_match": (
                len(g["indices"]) == nnz
            ),
            "matrix_attributes": {
                k: repr(v)[:300]
                for k, v in g.attrs.items()
            },
        })

        pointers_valid = (
            len(indptr) == n_cells + 1
            and indptr[0] == 0
            and indptr[-1] == nnz
            and np.all(indptr[:-1] <= indptr[1:])
        )

        result["sparse_pointers_valid"] = bool(pointers_valid)

        print("\n=== MATRIX INFORMATION ===")
        print("Genes:", n_genes)
        print("Nuclei:", n_cells)
        print("Stored entries:", nnz)
        print("Matrix shape:", (n_genes, n_cells))
        print("Data dtype:", g["data"].dtype)

        print("\nFirst genes:", genes[:10])
        print("First barcodes:", barcodes[:10])

        print("\nDuplicate genes:", result["gene_duplicates"])
        print("Duplicate barcodes:", result["barcode_duplicates"])
        print("Pointers valid:", pointers_valid)

        if not pointers_valid:
            result["error"] = "Invalid sparse matrix pointers"
            return result, genes

        result["stored_entries_per_nucleus"] = stats(
            np.diff(indptr)
        )

        # Inspect real expression values of first nucleus
        if n_cells > 0:
            start = int(indptr[0])
            end = min(int(indptr[1]), start + 10)

            gene_indices = g["indices"][start:end]
            counts = g["data"][start:end]

            examples = []

            print("\n=== FIRST NUCLEUS EXPRESSION ===")

            for idx, value in zip(gene_indices, counts):
                idx = int(idx)
                if 0 <= idx < n_genes:
                    examples.append({
                        "gene": genes[idx],
                        "count": int(value),
                    })
                    print(genes[idx], int(value))

            result["first_nucleus_expression"] = examples

        if not deep:
            return result, genes

        # Deep QC: process the sparse matrix in chunks
        print("\nRunning deep nucleus-level QC...")

        mitochondrial = np.asarray([
            gene.upper().startswith("MT-")
            for gene in genes
        ])

        qc_rows = []
        negative_values = 0
        stored_zeros = 0
        invalid_indices = 0

        for start_cell in range(0, n_cells, chunk_cells):

            end_cell = min(start_cell + chunk_cells, n_cells)

            start_idx = int(indptr[start_cell])
            end_idx = int(indptr[end_cell])

            values = g["data"][start_idx:end_idx]
            indices = g["indices"][start_idx:end_idx]

            per_cell_entries = np.diff(
                indptr[start_cell:end_cell + 1]
            ).astype(np.int64)

            column_ids = np.repeat(
                np.arange(end_cell - start_cell),
                per_cell_entries,
            )

            total_counts = np.bincount(
                column_ids,
                weights=values,
                minlength=end_cell - start_cell,
            )

            detected_genes = np.bincount(
                column_ids[values != 0],
                minlength=end_cell - start_cell,
            )

            valid = (indices >= 0) & (indices < n_genes)

            mt_mask = np.zeros(len(indices), dtype=bool)
            mt_mask[valid] = mitochondrial[indices[valid]]

            mt_counts = np.bincount(
                column_ids[mt_mask],
                weights=values[mt_mask],
                minlength=end_cell - start_cell,
            )

            negative_values += int(np.count_nonzero(values < 0))
            stored_zeros += int(np.count_nonzero(values == 0))
            invalid_indices += int(np.count_nonzero(~valid))

            for j in range(end_cell - start_cell):

                total = float(total_counts[j])

                pct_mt = (
                    100 * float(mt_counts[j]) / total
                    if total > 0 else 0.0
                )

                qc_rows.append((
                    barcodes[start_cell + j],
                    int(total),
                    int(detected_genes[j]),
                    round(pct_mt, 4),
                ))

        qc_file = output / (
            path.stem + ".nucleus_qc.csv"
        )

        with qc_file.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "barcode",
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_MT",
            ])
            writer.writerows(qc_rows)

        qc_array = np.array(
            [[r[1], r[2], r[3]] for r in qc_rows],
            dtype=float,
        )

        result["deep_qc"] = {
            "total_counts": stats(qc_array[:, 0]),
            "detected_genes": stats(qc_array[:, 1]),
            "pct_mitochondrial": stats(qc_array[:, 2]),
            "empty_nuclei": int(
                np.count_nonzero(qc_array[:, 0] == 0)
            ),
            "negative_values": negative_values,
            "stored_zero_values": stored_zeros,
            "invalid_gene_indices": invalid_indices,
            "qc_file": str(qc_file),
        }

        print("QC saved:", qc_file)

    return result, genes


def compare_huuki(path, nemo_genes):

    if not path.exists():
        return {
            "status": "Huuki file not found",
            "path": str(path),
        }

    try:
        import anndata as ad

        adata = ad.read_h5ad(path, backed="r")

        try:
            huuki_genes = set(map(str, adata.var_names))
            nemo_set = set(nemo_genes)

            common = nemo_set & huuki_genes

            comparison = {
                "huuki_nuclei": int(adata.n_obs),
                "huuki_genes": int(adata.n_vars),
                "nemo_genes": len(nemo_set),
                "shared_genes": len(common),
                "nemo_overlap_percent": round(
                    100 * len(common) / len(nemo_set), 2
                ),
                "huuki_overlap_percent": round(
                    100 * len(common) / len(huuki_genes), 2
                ),
                "huuki_obs_columns": list(adata.obs.columns),
                "huuki_var_columns": list(adata.var.columns),
                "note": (
                    "Gene matching uses var_names. "
                    "Check gene ID systems if overlap is low."
                ),
            }

        finally:
            adata.file.close()

        return comparison

    except Exception as e:
        return {
            "status": "Huuki comparison failed",
            "error": str(e),
        }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--max-files", type=int, default=0,
        help="Maximum donor files; 0 means all",
    )
    parser.add_argument(
        "--deep", action="store_true",
        help="Compute complete nucleus-level QC",
    )
    parser.add_argument(
        "--chunk-cells", type=int, default=256,
    )

    args = parser.parse_args()

    if args.chunk_cells <= 0:
        parser.error("--chunk-cells must be positive")

    OUTPUT.mkdir(parents=True, exist_ok=True)

    all_files = sorted(INPUT.glob("*.h5"))

    donor_files = [
        p for p in all_files
        if p.name.endswith(".gmg.donors.h5")
    ]

    categories = Counter(
        "donor_counts"
        if p.name.endswith(".gmg.donors.h5")
        else "metacells"
        if "metacells" in p.name
        else "digital_expression"
        if "digital_expression" in p.name
        else "other"
        for p in all_files
    )

    print("\n=== DATASET INVENTORY ===")
    print("All H5 files:", len(all_files))
    print("File categories:", dict(categories))
    print("Donor H5 files:", len(donor_files))

    if args.max_files > 0:
        donor_files = donor_files[:args.max_files]

    if not donor_files:
        print("No matching donor H5 files found.")
        return

    results = []

    reference_genes = None
    reference_hash = None

    for i, path in enumerate(donor_files, 1):

        print(f"\n[{i}/{len(donor_files)}] {path.name}")

        result, genes = inspect_file(
            path, OUTPUT,
            deep=args.deep,
            chunk_cells=args.chunk_cells,
        )

        if genes is not None:

            if reference_genes is None:
                reference_genes = set(genes)
                reference_hash = result["gene_order_sha256"]

                (OUTPUT / "first_batch_genes.txt").write_text(
                    "\n".join(genes) + "\n"
                )

            result["same_gene_order_as_first_batch"] = (
                result["gene_order_sha256"] == reference_hash
            )

            result["gene_overlap_with_first_batch"] = len(
                reference_genes & set(genes)
            )

        report_path = OUTPUT / (
            path.stem + ".inspection.json"
        )

        report_path.write_text(
            json.dumps(result, indent=2)
        )

        results.append(result)

    # Batch summary
    fields = [
        "filename",
        "size_MB",
        "genes",
        "nuclei",
        "stored_entries",
        "data_dtype",
        "gene_duplicates",
        "barcode_duplicates",
        "matrix_shape_consistent",
        "data_indices_lengths_match",
        "sparse_pointers_valid",
        "same_gene_order_as_first_batch",
        "gene_overlap_with_first_batch",
        "error",
    ]

    with (OUTPUT / "donor_files_summary.csv").open(
        "w", newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh, fieldnames=fields,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(results)

    # All-file inventory
    with (OUTPUT / "all_h5_filenames.csv").open(
        "w", newline=""
    ) as fh:

        writer = csv.writer(fh)
        writer.writerow(["filename", "size_MB"])

        for p in all_files:
            writer.writerow([
                p.name,
                round(p.stat().st_size / 1024**2, 2),
            ])

    # Huuki comparison
    huuki_comparison = (
        compare_huuki(HUUKI, list(reference_genes))
        if reference_genes else None
    )

    overall = {
        "input_directory": str(INPUT),
        "all_h5_files": len(all_files),
        "categories": dict(categories),
        "donor_files_available": len([
            p for p in all_files
            if p.name.endswith(".gmg.donors.h5")
        ]),
        "donor_files_inspected": len(results),
        "total_nuclei_across_batches": sum(
            r.get("nuclei", 0) for r in results
        ),
        "same_gene_order_all_batches": all(
            r.get("same_gene_order_as_first_batch", False)
            for r in results
        ),
        "deep_qc": args.deep,
        "huuki_comparison": huuki_comparison,
    }

    (OUTPUT / "overall_report.json").write_text(
        json.dumps(overall, indent=2)
    )

    print("\n========== COMPLETE ==========")
    print("Reports directory:", OUTPUT)
    print("Donor summary: donor_files_summary.csv")
    print("Overall report: overall_report.json")
    print("Huuki comparison:")
    print(json.dumps(huuki_comparison, indent=2))


if __name__ == "__main__":
    main()
