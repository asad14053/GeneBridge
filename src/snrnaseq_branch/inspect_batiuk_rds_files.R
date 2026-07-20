project_root <- "/users/mjabin/projects/GeneBridge"

batiuk_dir <- file.path(project_root, "data/raw/Batiuk")

files <- c(
  "annotations_final.RDS",
  "snRNA-seq_Conos_object.RDS",
  "snRNA-seq_raw_countmatrices.RDS",
  "visium_spatial_transcriptomics_raw_countmatrices.RDS"
)

out_dir <- file.path(project_root, "outputs/batiuk_inspection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in files) {
  path <- file.path(batiuk_dir, f)

  cat("\n============================================================\n")
  cat("FILE:", f, "\n")
  cat("PATH:", path, "\n")
  cat("EXISTS:", file.exists(path), "\n")

  if (!file.exists(path)) next

  obj <- readRDS(path)

  cat("CLASS:\n")
  print(class(obj))

  cat("\nOBJECT SUMMARY:\n")
  print(obj)

  cat("\nNAMES:\n")
  print(names(obj))

  sink(file.path(out_dir, paste0(f, "_structure.txt")))
  cat("FILE:", f, "\n")
  cat("CLASS:\n")
  print(class(obj))
  cat("\nSUMMARY:\n")
  print(obj)
  cat("\nNAMES:\n")
  print(names(obj))
  cat("\nSTRUCTURE:\n")
  str(obj, max.level = 3)
  sink()
}

cat("\nDONE. Outputs saved to:\n")
cat(out_dir, "\n")