#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

# Canonical N24 biological metadata
INPUT_METADATA = (
    ROOT
    / "data/metadata/xenium_sample_metadata_24_common.csv"
)

OUTPUT_DIR = ROOT / "data/metadata"
OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_METADATA = (
    OUTPUT_DIR
    / "xenium_DE_metadata_24.csv"
)

OUTPUT_REFERENCE = (
    OUTPUT_DIR
    / "xenium_PI_run_metadata_24.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "xenium_DE_metadata_24_validation.txt"
)


# =============================================================================
# PI Xenium slide/run metadata
# =============================================================================

groups = [

    (
        "output-XETG00133__0021210",
        "2024-11-18_xenium-pnn-SKCCC",
        ["Br5400", "Br5622", "Br6437"],
    ),

    (
        "output-XETG00133__0021218",
        "2024-11-18_xenium-pnn-SKCCC",
        ["Br5314", "Br5588", "Br5746"],
    ),

    (
        "output-XETG00133__0033849",
        "2024-12-06_xenium-pnn-SKCCC",
        ["Br5639", "Br8433", "Br8772"],
    ),

    (
        "output-XETG00133__0033862",
        "2024-12-06_xenium-pnn-SKCCC",
        ["Br1139", "Br5973", "Br8667"],
    ),

    (
        "output-XETG00133__0033839",
        "2024-12-12_xenium-pnn-SKCCC",
        ["Br5436", "Br5931", "Br6496"],
    ),

    (
        "output-XETG00133__0033841",
        "2024-12-12_xenium-pnn-SKCCC",
        ["Br2421", "Br6032", "Br6389"],
    ),

    (
        "output-XETG00133__0034238",
        "2024-11-12_Xenium-pnn-SKCCC",
        ["Br2039", "Br2719", "Br6432"],
    ),

    (
        "output-XETG00133__0034252",
        "2024-11-12_Xenium-pnn-SKCCC",
        ["Br1113", "Br5373", "Br5590"],
    ),
]


# =============================================================================
# Build slide/run reference
# =============================================================================

rows = []

for slide_id, run_date, donors in groups:

    for brnum in donors:

        rows.append(
            {
                "BrNum": brnum,
                "slide_id": slide_id,
                "run_date": run_date,
            }
        )


reference = pd.DataFrame(rows)

assert reference.shape == (24, 3)
assert reference["BrNum"].nunique() == 24
assert reference["slide_id"].nunique() == 8
assert reference["run_date"].nunique() == 4
assert not reference.isna().any().any()

reference = (
    reference
    .sort_values("BrNum")
    .reset_index(drop=True)
)

reference.to_csv(
    OUTPUT_REFERENCE,
    index=False,
)


# =============================================================================
# Biological metadata
# =============================================================================

metadata = pd.read_csv(
    INPUT_METADATA
)

required = {
    "BrNum",
    "Dx",
    "Age",
    "Sex",
}

missing = (
    required
    - set(metadata.columns)
)

if missing:

    raise RuntimeError(
        "N24 donor metadata missing columns: "
        f"{sorted(missing)}"
    )


metadata = (
    metadata[
        [
            "BrNum",
            "Dx",
            "Age",
            "Sex",
        ]
    ]
    .copy()
)

metadata["BrNum"] = (
    metadata["BrNum"]
    .astype(str)
)


# =============================================================================
# N24 integrity
# =============================================================================

assert metadata.shape[0] == 24
assert metadata["BrNum"].nunique() == 24
assert not metadata.duplicated("BrNum").any()
assert "Br6432" in set(metadata["BrNum"])

assert metadata["Dx"].isin(
    ["SCZ", "NTC"]
).all()


reference_donors = set(
    reference["BrNum"]
)

metadata_donors = set(
    metadata["BrNum"]
)

assert metadata_donors == reference_donors, (
    "Biological metadata and PI run metadata "
    "contain different donors.\n"
    f"Missing biological: "
    f"{sorted(reference_donors - metadata_donors)}\n"
    f"Unexpected biological: "
    f"{sorted(metadata_donors - reference_donors)}"
)


# =============================================================================
# Merge
# =============================================================================

metadata = metadata.merge(
    reference,
    on="BrNum",
    how="left",
    validate="one_to_one",
)

metadata = metadata[
    [
        "BrNum",
        "Dx",
        "Age",
        "Sex",
        "slide_id",
        "run_date",
    ]
]

assert metadata.shape == (24, 6)
assert metadata["BrNum"].nunique() == 24
assert not metadata.duplicated("BrNum").any()
assert not metadata.isna().any().any()

assert "Br6432" in metadata["BrNum"].values


# =============================================================================
# Critical slide check
# =============================================================================

slide_34238 = set(

    metadata.loc[
        metadata["slide_id"]
        == "output-XETG00133__0034238",
        "BrNum",
    ]
)

assert slide_34238 == {
    "Br2039",
    "Br2719",
    "Br6432",
}


# =============================================================================
# Write
# =============================================================================

metadata = (
    metadata
    .sort_values("BrNum")
    .reset_index(drop=True)
)

metadata.to_csv(
    OUTPUT_METADATA,
    index=False,
)


# =============================================================================
# Validation report
# =============================================================================

report = f"""
XENIUM DE METADATA N24 VALIDATION
=================================

FINAL STATUS: PASS

Input biological metadata:
{INPUT_METADATA}

Canonical DE metadata:
{OUTPUT_METADATA}

PI run metadata:
{OUTPUT_REFERENCE}

Total donors:
{len(metadata)}

Unique BrNum:
{metadata['BrNum'].nunique()}

Br6432 included:
{'Br6432' in metadata['BrNum'].values}

Diagnosis counts:
{metadata['Dx'].value_counts().sort_index().to_string()}

Missing values:
{metadata.isna().sum().to_string()}

Slide counts:
{metadata['slide_id'].value_counts().sort_index().to_string()}

Run-date counts:
{metadata['run_date'].value_counts().sort_index().to_string()}

0034238 donors:
{sorted(slide_34238)}

Expected:
['Br2039', 'Br2719', 'Br6432']

Excluded donors:
NONE

FINAL STATUS: PASS
"""

OUTPUT_REPORT.write_text(
    report.strip() + "\n"
)

print(report)

print("\nCanonical N24 metadata:\n")

print(
    metadata.to_string(
        index=False
    )
)

print("\nWritten:")
print(OUTPUT_METADATA)
print(OUTPUT_REFERENCE)
print(OUTPUT_REPORT)
