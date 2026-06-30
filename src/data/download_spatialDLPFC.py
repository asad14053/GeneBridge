from pathlib import Path
import subprocess
import textwrap
import sys

# If this file is GeneBridge/scripts/download_spatialDLPFC.py,
# project_root becomes GeneBridge/
PROJECT_ROOT = Path(__file__).resolve().parents[1]

VISIUM_DIR = PROJECT_ROOT / "data" / "processed" / "visium"
SNRNA_DIR = PROJECT_ROOT / "data" / "processed" / "snrnaseq"
TMP_R_SCRIPT = PROJECT_ROOT / "data" / "download_spatialDLPFC.R"

VISIUM_DIR.mkdir(parents=True, exist_ok=True)
SNRNA_DIR.mkdir(parents=True, exist_ok=True)

r_code = f"""
project_root <- "{PROJECT_ROOT.as_posix()}"

if (!requireNamespace("BiocManager", quietly = TRUE)) {{
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}}

pkgs <- c(
    "spatialLIBD",
    "HDF5Array"
)

for (pkg in pkgs) {{
    if (!requireNamespace(pkg, quietly = TRUE)) {{
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }}
}}

library(spatialLIBD)
library(HDF5Array)

visium_dir <- file.path(project_root, "data", "processed", "visium")
snrna_dir  <- file.path(project_root, "data", "processed", "snrnaseq")

dir.create(visium_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(snrna_dir, recursive = TRUE, showWarnings = FALSE)

message("Downloading Visium SpatialExperiment...")
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

message("Visium object:")
print(spe)

saveRDS(
    spe,
    file = file.path(visium_dir, "spatialDLPFC_Visium_spe.rds")
)

message("Saved Visium RDS to:")
message(file.path(visium_dir, "spatialDLPFC_Visium_spe.rds"))

message("Downloading snRNA-seq zip...")
sce_path_zip <- spatialLIBD::fetch_data(type = "spatialDLPFC_snRNAseq")

message("snRNA-seq zip path:")
message(sce_path_zip)

message("Unzipping snRNA-seq data...")
unzip(sce_path_zip, exdir = snrna_dir)

sce_h5_path <- file.path(snrna_dir, "sce_DLPFC_annotated")

message("Loading snRNA-seq HDF5-backed object...")
sce <- HDF5Array::loadHDF5SummarizedExperiment(sce_h5_path)

message("snRNA-seq object:")
print(sce)

saveRDS(
    sce,
    file = file.path(snrna_dir, "spatialDLPFC_snRNAseq_sce_h5backed.rds")
)

message("Saved snRNA-seq RDS pointer to:")
message(file.path(snrna_dir, "spatialDLPFC_snRNAseq_sce_h5backed.rds"))

message("DONE.")
"""

TMP_R_SCRIPT.write_text(textwrap.dedent(r_code))

print(f"Project root: {PROJECT_ROOT}")
print(f"Temporary R script written to: {TMP_R_SCRIPT}")

try:
    subprocess.run(
        ["Rscript", str(TMP_R_SCRIPT)],
        check=True
    )
except FileNotFoundError:
    print("ERROR: Rscript was not found. You need R installed and available in PATH.")
    sys.exit(1)
except subprocess.CalledProcessError as e:
    print("ERROR: R script failed.")
    print(e)
    sys.exit(e.returncode)

print("Download completed successfully.")