# GeneBridge
**Visium-guided Xenium gene imputation framework and dataset for reconstructing cell-resolved spatial transcriptomes in schizophrenia brain tissue**

A computational framework for integrating Visium, Xenium, PsychENCODE, and Perturb-seq datasets to investigate schizophrenia-associated molecular signatures in the human dorsolateral prefrontal cortex (dlPFC).

---

## Project Goals

- Integrate Visium and Xenium spatial transcriptomics datasets
- Develop Visium-guided Xenium gene expression imputation
- Identify schizophrenia-associated spatial molecular signatures
- Characterize layer-specific and cell-type-specific alterations
- Validate findings using PsychENCODE and Perturb-seq datasets

---

## Repository Structure

```text
SCZ-Visium-Xenium-Imputation/
│
├── README.md
├── LICENSE
├── .gitignore
├── environment.yml
├── requirements.txt
│
├── data/
│   ├── raw/
│   │   ├── visium/
│   │   ├── xenium/
│   │   └── psychencode/
│   │
│   ├── processed/
│   │   ├── visium/
│   │   ├── xenium/
│   │   └── integrated/
│   │
│   └── metadata/
│       ├── sample_metadata.csv
│       ├── donor_metadata.csv
│       └── feature_metadata.csv
│
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_qc_analysis.ipynb
│   ├── 03_shared_gene_analysis.ipynb
│   ├── 04_spatial_visualization.ipynb
│   ├── 05_imputation_benchmark.ipynb
│   └── 06_scz_analysis.ipynb
│
├── src/
│   │
│   ├── preprocessing/
│   │   ├── load_visium.py
│   │   ├── load_xenium.py
│   │   ├── qc_visium.py
│   │   ├── qc_xenium.py
│   │   └── normalization.py
│   │
│   ├── integration/
│   │   ├── shared_genes.py
│   │   ├── celltype_mapping.py
│   │   ├── spatial_alignment.py
│   │   └── pseudospot_generation.py
│   │
│   ├── imputation/
│   │   ├── baseline.py
│   │   ├── tangram_model.py
│   │   ├── graph_model.py
│   │   └── evaluation.py
│   │
│   ├── analysis/
│   │   ├── deg_analysis.py
│   │   ├── layer_analysis.py
│   │   ├── celltype_analysis.py
│   │   └── scz_signature_analysis.py
│   │
│   └── utils/
│       ├── config.py
│       ├── logger.py
│       └── plotting.py
│
├── configs/
│   ├── visium_config.yaml
│   ├── xenium_config.yaml
│   └── training_config.yaml
│
├── outputs/
│   ├── qc_reports/
│   ├── figures/
│   ├── models/
│   └── tables/
│
├── docs/
│   ├── project_plan.md
│   ├── data_dictionary.md
│   ├── methodology.md
│   └── experiment_log.md
│
└── tests/
    ├── test_visium_loader.py
    ├── test_xenium_loader.py
    └── test_imputation.py
```

---

## Workflow

### Phase 1 — Data Understanding

- Dataset inventory
- Metadata generation
- Sample tracking
- Shared gene identification

### Phase 2 — Quality Control

- Visium QC
- Xenium QC
- Filtering
- Normalization

### Phase 3 — Spatial Integration

- Shared gene mapping
- Cell-type annotation
- Spatial alignment
- Pseudo-spot generation

### Phase 4 — Baseline Imputation

- Ridge Regression
- Random Forest
- Holdout-gene prediction

### Phase 5 — Advanced Imputation

- Tangram
- Graph Neural Networks
- Spatial Transformer Models

### Phase 6 — Schizophrenia Analysis

- Differential expression
- Layer-specific analysis
- Cell-type analysis
- SCZ signature discovery

### Phase 7 — Validation

- PsychENCODE validation
- Perturb-seq validation
- External benchmarking

---

## Day 1 Checklist

- [ ] Read study and dataset documentation
- [ ] Build sample metadata table
- [ ] Inventory Visium samples
- [ ] Inventory Xenium samples
- [ ] Generate dataset summary
- [ ] Identify shared genes
- [ ] Create QC pipeline

---

## Expected Outputs

```text
outputs/
├── qc_reports/
├── figures/
├── models/
└── tables/
```

Example outputs:

- visium_qc_report.csv
- xenium_qc_report.csv
- shared_genes.csv
- imputation_metrics.csv
- scz_deg_results.csv

---

## Future Directions

- Cell-level transcriptome reconstruction
- Multi-modal spatial integration
- Schizophrenia biomarker discovery
- AI-powered spatial genomics agent

---

## Author

**Md Asaduzzaman Jabin**  
Postdoctoral Researcher  
Emory University
