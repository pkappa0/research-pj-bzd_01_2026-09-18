# C1 Balance / Ataxia: consistency of PLIF-clinical correspondence

既存C1 direction consistency summaryの表示のみを行ったstandalone figure。相関、中央値、min/maxの再計算は行わず、既存値をそのまま描画する。frozen PLIF/Clinical fingerprintは読み替え・変更しない。用途は記述的QCと仮説生成に限る。

## Input files

Source run: `runs/20260920_wobbling_related_phenotypes/`, repository commit `5c5cf06fa231b456f52448b0ed61f550e4ad9bb9`.

- `tables/c1_direction_consistency_summary.csv` → `data/c1_direction_consistency_summary.csv` (byte-identical snapshot).
- `tables/feature_display_key.csv` → `data/feature_display_key.csv` (byte-identical snapshot).
- Repository `ANALYSIS_CONTRACT.md` → `data/ANALYSIS_CONTRACT_source.md`.

`input_manifest.json` records exact source paths and SHA256. Prior runs are preserved. The original summary already contains per-feature sign counts, median, minimum and maximum across four prespecified C1 terms: Ataxia, Balance disorder, Coordination abnormal, Gait disturbance. This run does not read clinical values to calculate or optimize a correlation, cutoff, or score.

## Figure construction

- Horizontal forest-style plot, F01–F48 in original frozen order. No magnitude sorting or feature filtering.
- Point = existing `median_rho_across_c1_terms`; line = existing `min_rho` to `max_rho`. No confidence interval calculation.
- Spearman rho axis fixed at −1 to +1; dashed reference at zero.
- Horizontal divider between F21 and F22. F01–F21 = alpha1; F22–F48 = alpha2. No pooling or averaging across receptor blocks.
- Compact labels map original features to e.g. `F01 | α1 | γ2 ASP56`. `BZD_GAMMA2` is displayed as `γ2`, `BZD_SITE` as `site`; source residue numbering is preserved (leading zeroes removed only for display). These are structural features, not an assertion of causal or mechanistic residues.
- One uniform point/line color is used; there is no significance encoding. The optional sign-count strip is omitted to preserve readability.

### Right annotations

All annotations use the **existing count columns**, not new correlations:

- With all four terms estimable: positive count 4 → `4/4 +`; negative count 4 → `4/4 -`; positive count 3 → `3/4 +`; negative count 3 → `3/4 -`; otherwise `mixed`.
- With fewer than four estimable terms: `n/n estimable (of 4)` explicitly marks incomplete coverage. F29 and F48 have `2/2 estimable (of 4)` in this input. This does not imply that the other two terms are zero. The numerator/denominator refer to available estimable terms; the final parenthetical preserves the original four-term total.
- Undefined summary values remain NA, with no point/interval plotted. Zero rho is neither positive nor negative. Median/min/max are never replaced by zero or reconstructed from counts.

Annotations are descriptive direction labels, not significance thresholds, rankings, biomarkers or hit calls. Terms are related phenotypes, not independent biological replicates. The observational unit of the underlying correlations is drug. Intervals describe variation among clinical terms, not statistical uncertainty.

The figure includes the requested title, subtitle and footnote:

> C1 Balance / Ataxia: consistency of PLIF–clinical correspondence
>
> Median and range of drug-level Spearman rho across four prespecified C1 terms
>
> Observational unit = drug. Intervals show variation across clinical terms, not confidence intervals. Association/correspondence only; no causal or mechanistic inference.

## Reproduce / outputs

Requires Python, numpy, pandas and matplotlib. From repository root:

```bash
python runs/20260920_c1_direction_consistency/scripts/plot_c1_direction_consistency.py
```

The script can also use the bundled snapshots without network access. It checks source hashes, the 48-row order and label mapping, receptor partition, estimable/count consistency, valid rho range and min ≤ median ≤ max. No docking, PLIF, clinical ROR, Spearman or C1 summary recomputation is performed.

Outputs:

- `figures/c1_plif_direction_consistency.png` — 3240×5040 pixels, 240 dpi.
- `figures/c1_plif_direction_consistency.svg` — vector points/intervals and preserved text.
- `data/figure_display_annotations.csv` — original displayed numeric values plus compact labels and right annotations, as display metadata only.
- `QC_REPORT.json`, `input_manifest.json`.

The main figure and numerical equality with the archived summary are checked before delivery. No source tables or historical figures are overwritten.
