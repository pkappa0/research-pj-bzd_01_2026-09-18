Run ID: 20260917_1925_residue_class_heatmap
Parent Run: 20260917_1756_grouped_heatmap
Analysis Phase: residue-class-first interaction fingerprint visualization
Started: 2026-09-17T19:25:53+09:00
Completed: 2026-09-17T19:25:54+09:00

# RESIDUE_CLASS_GROUPED_HEATMAP_REPORT

This run adds a residue-class-first visualization. The previous receptor-first grouped heatmaps and the original frequency heatmap remain unchanged.

## Row hierarchy

Rows are ordered as residue class → residue/site → receptor → interaction type. Alpha1 and alpha2 are therefore adjacent for the same common-position site; there is no alpha1 block or alpha2 block. Class boundaries use thick separators and site boundaries use medium separators.

Observed classes in the parent matrix: Aromatic, Hydrophilic / polar, Charged.
Classes with no observed candidate feature rows: Hydrophobic / aliphatic, Other. Empty classes were not populated with artificial rows.

## Two heatmaps

The raw heatmap displays the unchanged 0–1 interaction frequency. The row-normalized heatmap displays each feature's raw frequency minus its four-drug mean. The normalized view is for relative pattern reading only; the raw matrix remains the numeric source.

## Labels

Labels use the compact form `gamma2 TYR58 | α1 | pi-stack` or `alpha HIS102 | α2 | hydrophobic`. Residue class is shown as a heading rather than repeated in every row.

## QC

Input feature rows: 25; raw frequency columns: diazepam, alprazolam, triazolam, zolpidem; raw values copied without transformation into results/tables/residue_class_grouped_frequency_matrix.csv.

## Files

results/figures/residue_class_grouped_frequency_heatmap.png; results/figures/residue_class_grouped_row_normalized_heatmap.png; results/tables/residue_class_grouped_frequency_matrix.csv; results/tables/residue_class_grouped_row_normalized_matrix.csv
