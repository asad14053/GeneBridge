#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import gzip

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

CURRENT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "envi_layer_pseudobulk"
    / "ENVI_measured300_donor_layer_pseudobulk.csv.gz"
)

META = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "envi_layer_pseudobulk"
    / "ENVI_donor_layer_pseudobulk_metadata.csv"
)

SEARCH_ROOTS = [
    ROOT / "data",
    ROOT / "outputs",
]


def section(x):
    print("\n" + "=" * 110)
    print(x)
    print("=" * 110)


section("CURRENT ENVI 300-GENE PSEUDOBULK")

print("File:", CURRENT)
print("Exists:", CURRENT.exists())

if CURRENT.exists():
    print("Size:", CURRENT.stat().st_size, "bytes")

    try:
        x = pd.read_csv(CURRENT, nrows=5)
        print("Columns:", x.columns[:15].tolist())
        print("Number of columns:", len(x.columns))
        print("\nFirst rows:")
        print(x.iloc[:, :8].to_string(index=False))
    except Exception as e:
        print("READ ERROR:", repr(e))


section("CURRENT ENVI PSEUDOBULK METADATA")

print("File:", META)
print("Exists:", META.exists())

if META.exists():
    m = pd.read_csv(META)
    print("Shape:", m.shape)
    print("Columns:", m.columns.tolist())
    print(m.head().to_string(index=False))


section("SEARCH FOR POSSIBLE CELL-LEVEL / GENE-LEVEL OOF ENVI FILES")

keywords = (
    "oof",
    "heldout",
    "held_out",
    "imput",
    "predict",
    "prediction",
    "envi",
)

extensions = {
    ".h5ad",
    ".csv",
    ".gz",
    ".tsv",
    ".parquet",
    ".npy",
    ".npz",
    ".pkl",
}


hits = []

for root in SEARCH_ROOTS:

    if not root.exists():
        continue

    for p in root.rglob("*"):

        if not p.is_file():
            continue

        name = p.name.lower()

        if not any(k in name for k in keywords):
            continue

        # Include common data formats.
        if (
            p.suffix.lower() not in extensions
            and not name.endswith(".csv.gz")
            and not name.endswith(".tsv.gz")
        ):
            continue

        try:
            size = p.stat().st_size
        except Exception:
            size = None

        hits.append(
            {
                "path": str(p.relative_to(ROOT)),
                "size_bytes": size,
            }
        )


hits_df = pd.DataFrame(hits)

if hits_df.empty:
    print("No candidate files found.")
else:
    hits_df = hits_df.sort_values(
        ["path"]
    )

    print("Candidate files:", len(hits_df))

    print(
        hits_df.to_string(
            index=False,
            max_rows=300
        )
    )


section("LIKELY OOF CANDIDATES")

if hits_df.empty:
    print("NONE")

else:

    candidate_mask = hits_df["path"].str.lower().str.contains(
        r"oof|heldout|held_out|prediction|predicted|imputed",
        regex=True,
    )

    candidates = hits_df[candidate_mask].copy()

    if candidates.empty:
        print("No filename explicitly contains OOF / heldout / prediction / imputed.")

    else:
        print(
            candidates.to_string(
                index=False,
                max_rows=300
            )
        )


section("CURRENT FILE-NAME INTERPRETATION")

print(
    """
The currently used file is named:

    ENVI_measured300_donor_layer_pseudobulk.csv.gz

The word 'measured300' strongly suggests that this matrix may contain the
original measured Xenium 300 genes rather than OOF predictions.

The previous comparison found:

    Original vs current ENVI:
        exact entries = 48,300 / 48,300
        Pearson       = 1
        Spearman      = 1
        MAE           = 0
        RMSE          = 0

Therefore this file must NOT be treated as an OOF-imputation result until its
upstream source is verified.
"""
)

print("FINAL STATUS: SOURCE INSPECTION COMPLETE")
