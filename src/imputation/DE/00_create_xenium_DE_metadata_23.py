#!/usr/bin/env python

from pathlib import Path
import pandas as pd

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

INPUT_DE_METADATA = (
    ROOT
    / "outputs/imputation_full/DE/before_imputation/"
      "original_xenium_pseudobulk_donor_metadata.csv"
)

OUTPUT_DIR = ROOT / "data/metadata"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_METADATA = OUTPUT_DIR / "xenium_DE_metadata_23.csv"
OUTPUT_REFERENCE = OUTPUT_DIR / "xenium_PI_run_metadata_24.csv"
OUTPUT_REPORT = OUTPUT_DIR / "xenium_DE_metadata_23_validation.txt"

groups = [
    ("output-XETG00133__0021210", "2024-11-18_xenium-pnn-SKCCC",
     ["Br5400", "Br5622", "Br6437"]),

    ("output-XETG00133__0021218", "2024-11-18_xenium-pnn-SKCCC",
     ["Br5314", "Br5588", "Br5746"]),

    ("output-XETG00133__0033849", "2024-12-06_xenium-pnn-SKCCC",
     ["Br5639", "Br8433", "Br8772"]),

    ("output-XETG00133__0033862", "2024-12-06_xenium-pnn-SKCCC",
     ["Br1139", "Br5973", "Br8667"]),

    ("output-XETG00133__0033839", "2024-12-12_xenium-pnn-SKCCC",
     ["Br5436", "Br5931", "Br6496"]),

    ("output-XETG00133__0033841", "2024-12-12_xenium-pnn-SKCCC",
     ["Br2421", "Br6032", "Br6389"]),

    ("output-XETG00133__0034238", "2024-11-12_Xenium-pnn-SKCCC",
     ["Br2039", "Br2719", "Br6432"]),

    ("output-XETG00133__0034252", "2024-11-12_Xenium-pnn-SKCCC",
     ["Br1113", "Br5373", "Br5590"]),
]

rows = []

for slide_id, run_date, donors in groups:
    for brnum in donors:
        rows.append({
            "BrNum": brnum,
            "slide_id": slide_id,
            "run_date": run_date,
        })

reference = pd.DataFrame(rows)

assert reference.shape == (24, 3)
assert reference["BrNum"].nunique() == 24
assert reference["slide_id"].nunique() == 8
assert reference["run_date"].nunique() == 4
assert not reference.isna().any().any()

reference = reference.sort_values("BrNum").reset_index(drop=True)
reference.to_csv(OUTPUT_REFERENCE, index=False)

de = pd.read_csv(INPUT_DE_METADATA)

required = {"patient_id", "Dx", "Age", "Sex"}
missing = required - set(de.columns)

if missing:
    raise RuntimeError(
        f"Existing DE metadata missing columns: {sorted(missing)}"
    )

de = (
    de[["patient_id", "Dx", "Age", "Sex"]]
    .rename(columns={"patient_id": "BrNum"})
    .copy()
)

de["BrNum"] = de["BrNum"].astype(str)

assert len(de) == 23
assert de["BrNum"].nunique() == 23
assert "Br6432" not in set(de["BrNum"])

reference_donors = set(reference["BrNum"])
analysis_donors = set(de["BrNum"])

missing_from_analysis = reference_donors - analysis_donors
unexpected_analysis = analysis_donors - reference_donors

assert missing_from_analysis == {"Br6432"}, (
    f"Expected only Br6432 excluded; found {sorted(missing_from_analysis)}"
)

assert not unexpected_analysis, (
    f"Unexpected donors: {sorted(unexpected_analysis)}"
)

metadata = de.merge(
    reference,
    on="BrNum",
    how="left",
    validate="one_to_one",
)

metadata = metadata[
    ["BrNum", "Dx", "Age", "Sex", "slide_id", "run_date"]
]

assert metadata.shape == (23, 6)
assert metadata["BrNum"].nunique() == 23
assert not metadata.duplicated("BrNum").any()
assert not metadata.isna().any().any()
assert metadata["Dx"].isin(["SCZ", "NTC"]).all()
assert "Br6432" not in metadata["BrNum"].values

slide_34238 = set(
    metadata.loc[
        metadata["slide_id"] == "output-XETG00133__0034238",
        "BrNum",
    ]
)

assert slide_34238 == {"Br2039", "Br2719"}

metadata = metadata.sort_values("BrNum").reset_index(drop=True)
metadata.to_csv(OUTPUT_METADATA, index=False)

report = f"""
XENIUM DE METADATA VALIDATION
=============================

FINAL STATUS: PASS

Input:
{INPUT_DE_METADATA}

Canonical metadata:
{OUTPUT_METADATA}

PI reference:
{OUTPUT_REFERENCE}

PI reference donors: {len(reference)}
Analysis donors: {len(metadata)}
Unique BrNum: {metadata['BrNum'].nunique()}

Excluded donor:
Br6432

Diagnosis counts:
{metadata['Dx'].value_counts().sort_index().to_string()}

Missing values:
{metadata.isna().sum().to_string()}

Slide counts:
{metadata['slide_id'].value_counts().sort_index().to_string()}

Run-date counts:
{metadata['run_date'].value_counts().sort_index().to_string()}

0034238 remaining donors:
{sorted(slide_34238)}

Expected:
['Br2039', 'Br2719']

FINAL STATUS: PASS
"""

OUTPUT_REPORT.write_text(report.strip() + "\n")

print(report)
print("\nCanonical metadata:\n")
print(metadata.to_string(index=False))

print("\nWritten:")
print(OUTPUT_METADATA)
print(OUTPUT_REFERENCE)
print(OUTPUT_REPORT)
