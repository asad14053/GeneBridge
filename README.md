# GeneBridge

**Reconstructing the cell-resolved spatial transcriptome of the human dlPFC in schizophrenia by integrating Xenium, Visium and snRNA-seq.**

> **Status: active research code, not a release.** The analyses are still being developed and scripts change often. Paths are currently hard-coded for the lab cluster, so this repository cannot yet run end-to-end anywhere else. See [Roadmap](#roadmap).

---

## Scientific question

Xenium measures single cells in place but covers only a **300-gene panel**. snRNA-seq covers the **whole transcriptome** but loses spatial position. Visium is **transcriptome-wide and spatial** but each spot mixes several cells.

GeneBridge asks:

> Can we reconstruct transcriptome-wide expression in spatially resolved Xenium cells, and does the reconstruction preserve cortical-layer organization, cell identity and schizophrenia (SCZ)-associated signal?

The project tests three things:

1. **Technical validity.** How accurately are held-out genes recovered?
2. **Biological validity.** Are spatial patterns, layers and cell-type specificity preserved?
3. **Scientific utility.** Does imputation change, or add to, SCZ differential-expression (DE) conclusions?

---

## Data

| Modality | Source | Role in GeneBridge |
|---|---|---|
| Xenium (300-gene panel) | 24 dlPFC donors (SCZ + neurotypical controls, NTC), matched to Visium | Imputation **target**; ground truth for held-out genes |
| Visium | spatialDLPFC (LIBD), the 24 matched donors | Independent spatial validation of imputed genes; layer annotation |
| snRNA-seq | Huuki-Myers et al. (10 brains), Batiuk et al., Ling et al. | Imputation **reference** |

Raw and processed data are **not** stored in this repository. They live under `data/` and `outputs/` on the cluster; both are git-ignored. See [`docs/xenium_output_files.md`](docs/xenium_output_files.md) for the Xenium file formats.

---

## Benchmark design

- **Gene-wise 3-fold cross-validation** over the 300 Xenium genes: each fold has 200 observed (training) genes and 100 held-out genes. Folds are canonical and shared by every method (`src/imputation/apply_canonical_gene_folds_N24.py`).
- **Leakage control.** Held-out genes are removed from `.X`, `.layers` and `.raw` before any method sees the data.
- **Methods:** ENVI, VISTA, gimVI, Tangram, SpaGE, TransImp.
- **Out-of-fold (OOF) stitching.** The three folds are combined into one 300-gene prediction per donor, so each gene is predicted exactly once (`<method>/01_stitch_oof_predictions.py`).
- **Metrics** (`src/imputation/common/benchmark_evaluation_v2.py`): Spearman correlation (SCC), SSIM, RMSE, MAE, JSD, and the absolute error in Moran's I.
- **Reference experiments:**

  | ID | snRNA-seq reference | Target |
  |---|---|---|
  | `ex5` | Br8667 only (donor-matched) | Br8667 Xenium |
  | `ex5_1` | All 10 Huuki brains | All 24 Xenium donors |
  | `ex5_3` | 9 Huuki brains, target donor **excluded** | Br8667 / Br6432 |

---

## Pipeline

| Stage | What it does | Main code |
|---|---|---|
| 1. Acquisition | Download Xenium/Visium (24 matched donors) and spatialDLPFC | `src/data/` |
| 2. Xenium processing | Build per-donor objects, merge N24, QC, cell-type (BANKSY) and layer annotation | `src/xenium_branch/` |
| 2. Visium processing | Convert RDS to h5ad, matched-donor layer annotation | `src/visium_branch/` |
| 2. snRNA-seq references | Huuki, Batiuk and Ling: conversion, gene harmonization, SCZ/NTC splits | `src/snrnaseq_branch/`, `scripts/ling_harmonization/` |
| 3. Benchmark | Folds, per-method fold runs, OOF stitching, evaluation | `src/imputation/{ENVI,VISTA,gimVI,Tangram,SpaGE,TransImp,common}/` |
| 3. Benchmark figures | Metric distributions, spatial gene panels, UMAP, NMI/ARI | `src/imputation/visualization/` |
| 4. Full transcriptome | ENVI production run on all 24 donors | `src/imputation/ENVI/run_envi_production.py`, `src/imputation/production/` |
| 5. Spatial validation | Cross-technology (Xenium vs Visium) comparison, Moran's I, three-source validation | `src/analysis/` |
| 5. Disease / DE | Layer- and cell-type pseudobulk, limma DE, concordance with published SCZ results | `src/imputation/DE/` |
| (future) AIVC | Data audits for a virtual-cell model | `src/AIVC/` |

[`docs/PIPELINE.md`](docs/PIPELINE.md) lists each SLURM job and the script it runs, stage by stage.

---

## Repository layout

```text
GeneBridge/
├── src/
│   ├── data/                 # downloads
│   ├── xenium_branch/        # Xenium build, QC, annotation
│   ├── visium_branch/        # Visium conversion and annotation
│   ├── snrnaseq_branch/      # snRNA-seq references (Huuki, Batiuk)
│   ├── imputation/
│   │   ├── <Method>/         # run_<method>_fold.py + 01_stitch_oof_predictions.py
│   │   ├── common/           # shared fold aggregation, evaluation, paired tests
│   │   ├── visualization/    # benchmark figures
│   │   ├── production/       # full-transcriptome targets
│   │   ├── reference_preparation/
│   │   └── DE/               # pseudobulk + differential expression
│   ├── analysis/             # cross-technology and spatial validation figures
│   ├── celltype_comparison/
│   └── AIVC/
├── scripts/
│   ├── slurm/                # cluster job scripts (imputation/, imputation/DE/, AIVC/)
│   ├── ling_harmonization/   # Ling et al. reference build
│   └── inspection/
├── docs/                     # pipeline map, VISTA patches, env snapshots
└── environment.yml           # minimal base environment (see below)
```

Files numbered `NN_` that are directly under `src/imputation/` (for example `03_align_vista_Br8667_genes.py`) come from the **Br8667 pilot**. The N24 production versions are in the method and `common/` folders.

---

## Environments

The jobs use four conda environments on the cluster:

| Env | Used for |
|---|---|
| `genebridge_scanpy` | Preprocessing, evaluation, figures, DE prep (Python) |
| `genebridge_r` | Visium/snRNA conversion, limma DE, enrichment (R) |
| `Imputation` | VISTA, gimVI, Tangram, SpaGE, TransImp (scvi-tools, PyTorch, torch-geometric) |
| `envi` | ENVI (`scenvi`) |

`environment.yml` currently only describes a minimal base env. Package snapshots of the `Imputation` env are in `software/` (kept on HGCC). Exported, pinned files for all four envs are on the roadmap.

---

## Third-party code

| Tool | How it is included | Notes |
|---|---|---|
| VISTA | Vendored in `src/imputation/VISTA/` | Locally patched; diffs in [`docs/patches/`](docs/patches/) |
| SpaGE | Vendored in `src/imputation/SpaGE/` | MIT license, © 2020 tabdelaal |
| TransImp (tranSpa) | Git submodule, `src/imputation/TransImp/tranSpa` | `git submodule update --init` |
| ENVI, Tangram, gimVI | Installed packages (`scenvi`, `tangram-sc`, `scvi-tools`) | — |

Please cite the original methods when using these wrappers.

---

## Running on the cluster

```bash
git clone --recurse-submodules https://github.com/asad14053/GeneBridge.git
cd GeneBridge
# data/ and outputs/ are expected as folders or symlinks at the repo root
sbatch scripts/slurm/imputation/run_envi_N24_fold_array.slurm   # example
```

Most scripts currently set `PROJECT_ROOT` to `/beegfs/labs/hulab/projects/mjabin/GeneBridge`. Change this before running elsewhere.

---

## Roadmap

- [x] 24-donor Xenium cohort (including the Br6432 restoration), layer and cell-type annotation
- [x] Six-method, 3-fold OOF benchmark on N24 (`benchmark_v.1`)
- [ ] Finish the Ling et al. snRNA-seq reference
- [ ] Final cross-method benchmark summary and figures
- [ ] Independent validation of full-transcriptome predictions (Visium, held-out donors)
- [ ] Measure how imputation changes SCZ DE conclusions
- [ ] Replace hard-coded paths with a single `config/paths.yaml`
- [ ] Pinned environment files for all four envs
- [ ] Reorganize into numbered pipeline stages and add a figure → script map for the manuscript

---

## Citation

Manuscript in preparation. Please contact the author before reusing the results.

## License

MIT. See [LICENSE](LICENSE). Vendored third-party code keeps its own license.