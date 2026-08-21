#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

INFILE = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "ENVI_10celltypes_imputed34687_DE_combined.csv.gz"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates"
)

OUT.mkdir(parents=True, exist_ok=True)

# NOTE:
# File currently has .gz extension but contains plain CSV text.
df = pd.read_csv(
    INFILE,
    compression=None
)

cols = [
    "gene",
    "P.value",
    "FDR_within_celltype",
    "FDR_global_10celltypes",
    "log2FC",
    "direction",
]


def extract(cell_type, expected_n=None):

    x = df[
        (df["cell_type"] == cell_type)
        &
        (df["FDR_within_celltype"] < 0.05)
    ].copy()

    x = x.sort_values(
        [
            "FDR_within_celltype",
            "P.value",
        ]
    ).reset_index(drop=True)

    x.insert(
        0,
        "rank",
        range(1, len(x) + 1)
    )

    if expected_n is not None:
        assert len(x) == expected_n, (
            f"{cell_type}: expected {expected_n}, "
            f"found {len(x)}"
        )

    print("\n" + "=" * 120)
    print(
        f"{cell_type} — ALL GENES WITH "
        f"WITHIN-CELL-TYPE FDR < 0.05"
    )
    print("=" * 120)

    print("Total         :", len(x))
    print(
        "Higher in SCZ :",
        int((x["log2FC"] > 0).sum())
    )
    print(
        "Higher in NTC :",
        int((x["log2FC"] < 0).sum())
    )

    print()

    print(
        x[
            ["rank"] + cols
        ].to_string(index=False)
    )

    slug = cell_type.replace("/", "_").replace(" ", "_")

    outfile = (
        OUT
        / f"ENVI_{slug}_withinFDR05_all.csv"
    )

    x[
        ["rank"] + cols
    ].to_csv(
        outfile,
        index=False
    )

    print("\nSaved:")
    print(outfile)

    return x


l6 = extract(
    "L6 Ex",
    expected_n=141
)

l5 = extract(
    "L5 Ex",
    expected_n=2
)


print("\n" + "=" * 120)
print("FINAL SUMMARY")
print("=" * 120)

print(
    f"L6 Ex: {len(l6)} total | "
    f"{(l6['log2FC'] > 0).sum()} higher SCZ | "
    f"{(l6['log2FC'] < 0).sum()} higher NTC"
)

print(
    f"L5 Ex: {len(l5)} total | "
    f"{(l5['log2FC'] > 0).sum()} higher SCZ | "
    f"{(l5['log2FC'] < 0).sum()} higher NTC"
)
