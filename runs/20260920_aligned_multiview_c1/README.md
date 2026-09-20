# Aligned multi-view structural–clinical figure

Visualization only, using existing frozen values and saved correlations. Historical runs and ANALYSIS_CONTRACT.md are unchanged; the contract is snapshotted here. Same receptor → different drugs remains the comparison unit. No docking, PLIF, FAERS ROR or correlation is recomputed. No composite score, feature selection, or clinical optimization is performed.

## Sources

All paths below are relative to repository root. Byte-identical snapshots are in this run; source_manifest.json provides SHA-256 hashes.

| Snapshot | Original source |
|---|---|
| `data/plif.csv` | `runs/20260920_vina_score_layer/data/frozen_plif_fingerprint.csv` |
| `data/clinical.csv` | `runs/20260920_vina_score_layer/data/frozen_clinical_fingerprint.csv` |
| `data/vina_summary.csv` | `runs/20260920_vina_score_layer/data/vina_scores_drug_receptor_summary.csv` |
| `data/plif_c1_correspondence.csv` | `runs/20260920_wobbling_focused_correspondence/tables/wobbling_plif_clinical_correspondence.csv` |
| `data/vina_c1_correspondence.csv` | `runs/20260920_vina_score_layer/data/vina_clinical_c1_correspondence.csv` |
| `data/feature_display_key.csv` | `runs/20260920_wobbling_focused_correspondence/tables/feature_display_key.csv` |
| `data/drug_display_order.csv` | `runs/20260920_vina_score_layer/config/drug_display_order.csv` |
| `ANALYSIS_CONTRACT_source.md` | `ANALYSIS_CONTRACT.md` |

## Display specification

Panel A rows follow the frozen order: diazepam, alprazolam, triazolam, zolpidem, lorazepam, clonazepam, midazolam, temazepam, zopiclone, zaleplon. All three heatmaps share identical row boundaries.

- A1 displays frozen contact frequency (0–1), F01–F48 without selection/reordering. F01–F21 are alpha1; F22–F48 are alpha2, separated by a line. Exact residue/source labels are in data/feature_display_key.csv.
- A2 colors encode **within-receptor descriptive favorability ranks only**. Rank the negative of the existing median seed-best Vina score separately across ten drugs for each receptor: 1 = least favorable, 10 = most favorable, average ranks for ties. No z-score or clinical-informed transformation is used. data/vina_display_ranks.csv saves this display-only transform. Cells show the original median Vina docking score to three decimals (kcal/mol); full-precision scores and seed counts remain in data/vina_summary.csv. All conditions have five independent seeds. Colors compare relative position within each receptor and do not imply absolute affinity equivalence between receptor blocks.
- A3 displays existing C1 logROR for Ataxia, Balance disorder, Coordination abnormal and Gait disturbance. Three missing cells remain gray, labeled NA; no imputation. The displayed sequential scale starts at zero and ends at the maximum of these existing C1 values.
- B1 selects the four prespecified C1 columns from the existing PLIF correspondence table, preserving all 48 features. B2 uses existing Vina favorability rho. Both share a fixed −1 to +1 Spearman color scale. Cell text is the saved n_effective, including where rho is non-estimable (gray). Values are only pivoted/reindexed for plotting, never recalculated. B1 and B2 deliberately have different row counts and cell heights; their rows do not represent matching units.

The only new numeric transformation is the explicitly requested Vina within-receptor display ranking. It does not replace scores or affect any correlation.

## Framing

PLIF interaction pattern and Vina docking score remain complementary, independent views aligned by drug identity. This figure supports descriptive inspection; it does not determine which representation is superior, establish complementary predictive performance, rank receptor dominance, or assign causal effects to residues or binding. No fitted predictor or model comparison is performed.

Vina docking score is a computational scoring proxy, not experimental affinity. All correspondence analyses use drug as the observational unit. Association/correspondence only; not causation or mediation. Source limitations remain: frozen all-role FAERS logROR is observational and subject to reporting confounding; the alpha2 construct remains provisional, with source composition/structural caveats in ANALYSIS_CONTRACT_source.md. Seeds are computational searches, not biological replicates.

## Reproduction and outputs

From repository root, using Python, numpy, pandas and matplotlib:

```sh
python runs/20260920_aligned_multiview_c1/scripts/plot_multiview.py
```

Outputs:

- figures/multiview_plif_vina_c1_profile.png (4840 × 4620 pixels, 220 dpi)
- figures/multiview_plif_vina_c1_profile.svg (editable text)
- data/vina_display_ranks.csv
- source_manifest.json and QC_REPORT.json

Assertions verify drug order, feature order, receptor boundaries, matrix shapes and source snapshot hashes before/after plotting. The rendered PNG was visually reviewed for row alignment, legend scales, labels, NA cells and annotation spacing. All historical tracked files remain unchanged.
