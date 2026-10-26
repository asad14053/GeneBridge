# Pipeline map: SLURM job → script

Generated from the `scripts/slurm/` files. It lists every cluster job and the script(s) it runs, grouped by pipeline stage.
Scripts marked ⚠️ are called by a job but are **not in the repository**.

Jobs whose names contain `ex5`, `ex5_1` or `ex5_3` are the reference experiments described in the main README. `N24` jobs run all 24 Xenium donors.

## 0. Br8667 pilot (historical)

| SLURM job | Runs |
|---|---|
| `align_vista_Br8667_genes.slurm` | `src/imputation/03_align_vista_Br8667_genes.py` |
| `audit_expression_layers_experiments_5.slurm` | `src/imputation/08_audit_expression_layers_experiments_5.py` |
| `audit_expression_layers_lowmem.slurm` | `src/imputation/08_audit_expression_layers_lowmem.py` |
| `check_three_reference_target_alignment.slurm` | `src/imputation/07_check_three_reference_target_alignment.py` |
| `create_experiment_5_3_huuki_9brains.slurm` | `src/imputation/05_3_create_huuki_9brains_excluding_xenium.py` |
| `create_experiment_5_3_huuki_9samples.slurm` | `src/imputation/05_3_create_huuki_9samples_excluding_xenium.py` ⚠️ |
| `extract_xenium_Br8667.slurm` | `src/imputation/04_extract_xenium_Br8667.py` |
| `inspect_xenium_Br8667.slurm` | `src/imputation/05_inspect_xenium_Br8667.py` |
| `make_huuki_visium_Br8667_h5ad.slurm` | `src/imputation/02_make_huuki_visium_Br8667_h5ad.py` |
| `prepare_vista_xenium_Br8667.slurm` | `src/imputation/06_prepare_vista_xenium_Br8667.py` |

## 1. Data acquisition

| SLURM job | Runs |
|---|---|
| `download_24_common_xenium_parallel.slurm` | `src/data/download_24_common_xenium_visium.py` |
| `download_24_common_xenium_visium_array.slurm` | `src/data/download_24_common_xenium_visium.py` |
| `download_spatialDLPFC.slurm` | `src/data/download_spatialDLPFC.R` |

## 2a. Xenium processing & annotation

| SLURM job | Runs |
|---|---|
| `build_xenium_N24_raw_old.slurm` | `src/xenium_branch/01_build_datastructure/01_build_xenium_N24_raw.py` |
| `get_cell_layer_annotation.slurm` | `src/xenium_branch/03_annotation/02_banksy_celltype_annotation.py` |
| `layerwise_annotation.slurm` | `src/xenium_branch/03_annotation/01_layerwise_annotation.py` |
| `merge_xenium_N24_raw.slurm` | `src/xenium_branch/01_build_datastructure/02_merge_xenium_N24_raw.py` |
| `xenium_N24_raw_array.slurm` | `src/xenium_branch/01_build_datastructure/01b_build_xenium_donor_raw.py` |
| `xenium_qc.slurm` | `src/xenium_branch/02_qc/01_xenium_qc.py` |

## 2b. Visium processing

| SLURM job | Runs |
|---|---|
| `convert_matched_visium_rds_to_h5ad.slurm` | `src/visium_branch/03_annotation/02_convert_matched_visium_rds_to_h5ad.R` |
| `convert_spatialDLPFC_light_RDS_to_h5ad.slurm` | `src/visium_branch/01_build_datastructure/02_convert_spatialDLPFC_light_RDS_to_h5ad.R` |
| `visualize_matched_N24_visium_layers.slurm` | `src/visium_branch/03_annotation/01_visualize_matched_N24_visium_layers.R` |

## 2c. snRNA-seq references

| SLURM job | Runs |
|---|---|
| `chunk_huuki_sce_to_h5ad_full_allgenes.slurm` | `src/snrnaseq_branch/chunk_huuki_sce_to_h5ad_full_allgenes.R` |
| `compare_batiuk_huuki_genes.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/02_compare_batiuk_huuki_genes.py` |
| `create_baituk_huuki_overlap_h5ads.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/04_create_overlap_filtered_h5ads.py` |
| `create_br8667_snrna.slurm` | `src/snrnaseq_branch/create_brain_snrna_from_huuki_h5ad.py` |
| `export_huuki_snrna_metadata.slurm` | `src/snrnaseq_branch/02_metadata_qc/01_export_huuki_snrna_metadata.R` |
| `extract_baituk_conos_annotations.slurm` | `src/snrnaseq_branch/04_baituk_celltype_annotation/01_extract_conos_annotations.R` |
| `extract_baituk_genes.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/01_extract_batiuk_genes.R` |
| `filter_baituk_from_rds.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/05_filter_baituk_from_rds.R` |
| `filter_baituk_h5ad.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/05_filter_baituk_h5ad.py` |
| `filter_huuki_h5ad.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/06_filter_huuki_h5ad.py` |
| `ling_build_gene_crosswalk.slurm` | `src/analysis/ling_harmonization/01_build_ling_h5_gene_crosswalk.py` |
| `ling_build_huuki_compatible_reference.slurm` | `scripts/ling_harmonization/05_build_ling_huuki_compatible_reference.py` |
| `ling_huuki_kwon_gene_harmonization.slurm` | `scripts/ling_harmonization/02_harmonize_ling_huuki_kwon_genes.py` |
| `ling_reference_01_prepare.slurm` | `scripts/ling_harmonization/06_prepare_ling_array_reference.py` |
| `ling_reference_02_array.slurm` | `scripts/ling_harmonization/07_process_ling_library.py` |
| `ling_reference_03_merge.slurm` | `scripts/ling_harmonization/08_merge_ling_reference.py` |
| `make_baituk_full_h5ad.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/03_make_baituk_full_h5ad.R` |
| `make_baituk_scz_ntc_h5ad.slurm` | `src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/03_make_baituk_scz_ntc_h5ad.R` |
| `split_baituk_rds_scz_ntc.slurm` | `src/snrnaseq_branch/concat_baituk_samples_scz_ntc.py`<br>`src/snrnaseq_branch/export_baituk_rds_samples_to_h5ad.R` |
| `split_baituk_scz_ntc.slurm` | `src/snrnaseq_branch/split_baituk_scz_ntc_chunked.py` |

## 3. Imputation benchmark (folds, runs, stitching, evaluation)

| SLURM job | Runs |
|---|---|
| `create_balanced_gene_folds_h5ad.slurm` | `src/imputation_beta/01_create_balanced_gene_folds_h5ad.py` |
| `create_balanced_gene_folds_h5ad_N24.slurm` | `src/imputation/apply_canonical_gene_folds_N24.py` |
| `filter_common_nonempty_cells_N24.slurm` | `src/imputation/filter_common_nonempty_cells_N24.py` |
| `imputation/aggregate_imputation_3fold_v2.slurm` | `src/imputation/common/aggregate_three_folds_v2.py`<br>`src/imputation/common/benchmark_evaluation_v2.py` |
| `imputation/aggregate_vista_ex5.slurm` | `src/imputation/common/aggregate_three_folds.py` |
| `imputation/aggregate_vista_ex5_1.slurm` | `src/imputation/common/aggregate_three_folds.py` |
| `imputation/aggregate_vista_ex5_3.slurm` | `src/imputation/common/aggregate_three_folds.py` |
| `imputation/analyze_gimvi_experiments.slurm` | `src/imputation/common/advisor_gene_panel_imputation.py`<br>`src/imputation/common/paired_tests_imputation_experiments.py` |
| `imputation/analyze_imputation_model.slurm` | `src/imputation/common/compare_advisor_genes_across_experiments.py`<br>`src/imputation/common/paired_tests_imputation_experiments.py` |
| `imputation/compare_advisor_genes_across_experiments.slurm` | `src/imputation/common/compare_advisor_genes_across_experiments.py` |
| `imputation/filter_common_nonempty_cells.slurm` | `src/imputation/common/filter_common_nonempty_cells.py` |
| `imputation/plot_oof_advisor_genes.slurm` | `src/imputation/common/benchmark_evaluation.py`<br>`src/imputation/common/plot_oof_advisor_genes.py` |
| `imputation/prepare_ex5_1.slurm` | `src/imputation/common/prepare_ex5_1_manifest.py` |
| `imputation/run_envi_3fold_array.slurm` | `src/imputation/ENVI/run_envi_fold.py` |
| `imputation/run_envi_N24_3fold_array.slurm` | `src/imputation/ENVI/run_envi_fold.py` |
| `imputation/run_envi_N24_donor_array.slurm` | `src/imputation/ENVI/run_envi_fold.py` |
| `imputation/run_envi_N24_fold_array.slurm` | `src/imputation/ENVI/run_envi_fold.py` |
| `imputation/run_envi_ex5_6task_array.slurm` | `src/imputation/ENVI/run_envi_fold.py` |
| `imputation/run_gimvi_3fold_array.slurm` | `src/imputation/gimVI/run_gimvi_fold.py` |
| `imputation/run_gimvi_N24_fold_array.slurm` | `src/imputation/gimVI/run_gimvi_fold.py` |
| `imputation/run_spage_3fold_array.slurm` | `src/imputation/SpaGE/SpaGE/main.py`<br>`src/imputation/SpaGE/run_spage_fold.py` |
| `imputation/run_spage_N24_fold_array.slurm` | `src/imputation/SpaGE/run_spage_fold.py` |
| `imputation/run_tangram_3fold_array.slurm` | `src/imputation/Tangram/run_tangram_fold.py` |
| `imputation/run_tangram_N24_fold_array.slurm` | `src/imputation/Tangram/run_tangram_fold.py` |
| `imputation/run_transimp_N24_fold_array.slurm` | `src/imputation/TransImp/run_transimpspa_fold.py` |
| `imputation/run_transimpspa_3fold_array.slurm` | `src/imputation/TransImp/run_transimpspa_fold.py` |
| `imputation/run_vista_N24_donor_array.slurm` | `src/imputation/VISTA/run_vista_fold.py` |
| `imputation/run_vista_N24_fold_array.slurm` | `src/imputation/VISTA/run_vista_fold.py` |
| `imputation/run_vista_ex5_1_3fold_array.slurm` | `src/imputation/VISTA/run_vista_fold.py` |
| `imputation/run_vista_ex5_3_3fold_array.slurm` | `src/imputation/VISTA/run_vista_fold.py` |
| `imputation/run_vista_ex5_3fold_array.slurm` | `src/imputation/VISTA/run_vista_fold.py` |
| `imputation/run_vista_ex5_6task_array.slurm` | `src/imputation/VISTA/run_vista_fold.py` |
| `imputation/run_vista_paired_tests.slurm` | `src/imputation/common/paired_tests_vista_experiments.py` |
| `imputation/stitch_envi_oof_all.slurm` | `src/imputation/ENVI/01_stitch_oof_predictions.py` |
| `imputation/stitch_gimvi_oof_all.slurm` | `src/imputation/gimVI/01_stitch_oof_predictions.py` |
| `imputation/stitch_spage_oof_all.slurm` | `src/imputation/SpaGE/01_stitch_oof_predictions.py` |
| `imputation/stitch_tangram_oof_all.slurm` | `src/imputation/Tangram/01_stitch_oof_predictions.py` |
| `imputation/stitch_transimp_oof_all.slurm` | `src/imputation/TransImp/01_stitch_oof_predictions.py` |
| `imputation/stitch_vista_oof_all.slurm` | `src/imputation/VISTA/01_stitch_oof_predictions.py` |

## 4. Benchmark figures

| SLURM job | Runs |
|---|---|
| `imputation/compare_methods_by_expression.slurm` | `src/imputation/visualization/05_compare_methods_by_expression.py` |
| `imputation/plot_all_models_umap.slurm` | `src/imputation/visualization/10_umap_method_across_experiments.py` |
| `imputation/plot_metric_distributions_across_experiments.slurm` | `src/imputation/visualization/03_metric_distributions_across_experiments.py` |
| `imputation/plot_metrics_by_expression.slurm` | `src/imputation/visualization/04_metric_distributions_by_expression.py` |
| `imputation/plot_rmse_mae_logscale.slurm` | `src/imputation/visualization/06_rmse_mae_logscale_comparison.py` |
| `imputation/plot_spatial_11genes_across_experiments.slurm` | `src/imputation/visualization/02_spatial_11genes_across_experiments.py` |
| `imputation/plot_spatial_11genes_across_methods.slurm` | `src/imputation/visualization/01_spatial_11genes_across_methods.py` |
| `imputation/plot_vista_umap.slurm` | `src/imputation/visualization/09_umap_vista_across_experiments.py` |

## 5. Full-transcriptome production

| SLURM job | Runs |
|---|---|
| `imputation/create_target_specific_snrna_refs.slurm` | `src/imputation/reference_preparation/create_target_specific_snrna_references.py` |
| `imputation/extract_ex2_scz_targets.slurm` | `src/imputation/production/01_extract_ex2_scz_targets.py` |
| `imputation/prepare_ex1_ex2_targets.slurm` | `src/imputation/production/00_prepare_ex1_ex2_targets.py` |
| `imputation/prepare_ex1_ex2_targets_fast.slurm` | `src/imputation/production/00_prepare_ex1_ex2_targets_fast.py` |
| `imputation/run_envi_production_24donors.slurm` | `src/imputation/ENVI/run_envi_production.py` |
| `imputation/run_envi_production_Br2039.slurm` | `src/imputation/ENVI/run_envi_production.py` |

## 6. Spatial & cross-technology validation

| SLURM job | Runs |
|---|---|
| `Br8667_OOF300_vs_34K_gene_variation.slurm` | `src/analysis/19_Br8667_OOF300_vs_34K_gene_variation.py` |
| `Br8667_advisor_genewise_OOF_figure.slurm` | `src/analysis/20_Br8667_advisor_genewise_OOF_figure.py` |
| `Br8667_advisor_moran_cosine_comparison.slurm` | `src/analysis/21_Br8667_advisor_moran_cosine_comparison.py` |
| `crosstech_measured300_all24.slurm` | `src/analysis/03_crosstech_measured300.py` |
| `envi_PI_figure_12NTC.slurm` | `src/analysis/08_envi_fulltranscriptome_PI_figure.py` |
| `envi_centered_cosine_figures_12NTC.slurm` | `src/analysis/08_envi_fulltranscriptome_PI_figure.py` |
| `envi_centered_cosine_metrics_12NTC.slurm` | `src/analysis/06_envi_fulltranscriptome_crosstech.py` |
| `envi_full_crosstech_12NTC.slurm` | `src/analysis/06_envi_fulltranscriptome_crosstech.py` |
| `envi_full_crosstech_test_Br8667.slurm` | `src/analysis/06_envi_fulltranscriptome_crosstech.py` |
| `envi_moran_12NTC.slurm` | `src/analysis/09_envi_moran_visium_vs_imputed.py` |
| `envi_moran_fullgenes_12NTC.slurm` | `src/analysis/09_envi_moran_visium_vs_imputed.py` |
| `measured300_moran_12NTC.slurm` | `src/analysis/11_measured300_moran_positive_control.py` |
| `measured300_moran_scatter_12NTC.slurm` | `src/analysis/13_plot_measured300_patterned_scatter.py` |
| `run_gene_compare.slurm` | `src/analysis/gene_compare_one_validation_gene.py` ⚠️ |
| `three_source_Br8667_all_figures.slurm` | `src/analysis/16_plot_Br8667_B_vs_C_template.py` ⚠️<br>`src/analysis/17_plot_Br8667_A_vs_B_A_vs_C_template.py` ⚠️ |
| `three_source_Br8667_all_template.slurm` | `src/analysis/18_plot_Br8667_three_source_all_template.py` |

## 7. Pseudobulk & differential expression

| SLURM job | Runs |
|---|---|
| `imputation/10a_build_envi_celltype_pseudobulk_23donors.slurm` | `src/imputation/DE/paper_concordance/10a_build_envi_celltype_pseudobulk_donor.py` |
| `imputation/11_celltype_global_paper_before_after_DE.slurm` | `src/imputation/DE/paper_concordance/11_celltype_global_paper_before_after_DE.R` |
| `imputation/DE/02a_envi_layer_pseudobulk_23donor_array.slurm` | `src/imputation/DE/paper_concordance/02a_build_envi_layer_pseudobulk_donor.py` |
| `imputation/DE/02b_aggregate_envi_layer_pseudobulk.slurm` | `src/imputation/DE/paper_concordance/02b_aggregate_envi_layer_pseudobulk.py` |
| `imputation/DE/03_layer_adjusted_de_array.slurm` | `src/imputation/DE/paper_concordance/03_layer_adjusted_de.R` |
| `imputation/DE/task2_3_build_oof_pseudobulk.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/01_build_model_oof_layer_pseudobulk.py` |
| `imputation/DE/task2_3_celltype_adjusted_de_all_models.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/05_celltype_adjusted_de_all_models.R` |
| `imputation/DE/task2_3_celltype_adjusted_de_all_models_nofilter.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/05_celltype_adjusted_de_all_models_nofilter.R` |
| `imputation/DE/task2_3_celltype_adjusted_de_all_models_nofilter_noslide.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/05_celltype_adjusted_de_all_models_nofilter_noslide.R` |
| `imputation/DE/task2_3_layer_adjusted_de_all300_filterQC.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/02_layer_adjusted_de_all_models_all300_filterQC.R` |
| `imputation/DE/task2_3_layer_adjusted_de_all_models.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/02_layer_adjusted_de_all_models.R` |
| `imputation/DE/task2_3_layer_adjusted_de_all_models_nofilter.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/02_layer_adjusted_de_all_models_nofilter.R` |
| `imputation/DE/task2_3_layer_adjusted_de_all_models_nofilter_noslide.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/02_layer_adjusted_de_all_models_nofilter_noslide.R` |
| `imputation/DE/task2_3_layer_adjusted_de_all_models_nofilter_rundate.slurm` | `src/imputation/DE/task2_3_ntc_dysregulation/02_layer_adjusted_de_all_models_nofilter_rundate.R` |
| `imputation/build_envi_celltype_pseudobulk.slurm` | `src/imputation/DE/05_build_envi_celltype_pseudobulk.py` |
| `imputation/build_envi_full_pseudobulk.slurm` | `src/imputation/DE/02_build_envi_full_pseudobulk.py` |
| `imputation/run_L6_pathway_enrichment.slurm` | `src/imputation/DE/09_L6_pathway_enrichment.R` |
| `imputation/run_envi_celltype_de.slurm` | `src/imputation/DE/06_envi_celltype_limma_de.R` |
| `imputation/run_envi_transcriptome_de.slurm` | `src/imputation/DE/03_envi_transcriptome_limma_de.R` |

## 8. Utilities

| SLURM job | Runs |
|---|---|
| `check_gpu_wait_50g.sh` | _(generic launcher)_ |
| `imputation/submit_envi_N24_all.sh` | _(generic launcher)_ |
| `imputation/submit_vista_ex5_3fold.sh` | _(generic launcher)_ |
| `run_notebook_cpu.slurm` | _(generic launcher)_ |
| `run_notebook_gpu200b.slurm` | _(generic launcher)_ |
| `run_notebook_gpu200b_ipynb.slurm` | _(generic launcher)_ |
| `run_notebook_gpu200b_test.slurm` | _(generic launcher)_ |
| `run_notebook_gpu_fast50g.slurm` | _(generic launcher)_ |
| `run_notebook_gpu_rtx.slurm` | _(generic launcher)_ |
| `run_python_cpu.slurm` | _(generic launcher)_ |
| `submit_split_baituk_scz_ntc.sh` | _(generic launcher)_ |

## 9. AIVC (exploratory)

| SLURM job | Runs |
|---|---|
| `AIVC/00_aggregate_presplit_audit.slurm` | `src/AIVC/data_prep/00_aggregate_presplit_audit.py` |
| `AIVC/00_presplit_audit.slurm` | `src/AIVC/data_prep/00_presplit_audit.py` |
| `AIVC/00_presplit_audit_array.slurm` | `src/AIVC/data_prep/00_presplit_audit_donor.py` |

## Not yet wired to a SLURM job

Scripts run interactively or not yet scheduled include `src/analysis/04,05,07,10,12,15_*`, most of `src/imputation/DE/paper_concordance/` (except 02a, 02b, 03, 10a, 11), `src/imputation/visualization/07,08,11_*`, and `src/imputation/validate_all_models_oof.py`.
