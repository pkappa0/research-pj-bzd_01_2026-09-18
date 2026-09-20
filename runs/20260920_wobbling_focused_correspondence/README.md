# Wobbling-focused PLIF-clinical phenotype correspondence

**association/correspondence only; not causation or mediation**

「ふらつき / motor instability」に焦点を当て、既存frozen X/Yと既存drug-level Spearman結果を表示した図。新しいendpointの合成、docking、PLIF、ROR、Spearmanの再計算やfeature selectionは行わない。統計単位は薬剤（最大n=10）。利用者指定の表示グループ化であり、元の14-term clinical fingerprintを変更しない。

## Panels and interpretation

- **A — Clinical phenotype definition**: PRIMARYはC1 Ataxia / Balance disorder / Coordination abnormal / Gait disturbance。SECONDARYはC5 Muscular weakness / Hypotoniaで、C1と同一視しない。C3 Somnolence / Sedation、C4 Dizziness / Vertigoは「contributing / potentially confounding phenotypes」、C2 Fallは「downstream outcome」とする文脈上の配置。これは概念整理であり因果順序を検証したものではない。矢印は使用しない。
- **B — Drug-level PLIF correspondence with balance / motor phenotypes**: 全48 featureをF01–F48の元の列順で保持。C1の4語、区切り、C5の2語の順。色は既存Spearman rho（−1〜+1、0中心）、数字は既存n_effective。灰色はNA/undefined correlationでゼロではない。行はstructural PLIF featuresであり「原因残基」ではない。
- **C — Representative structural patterns**: F02/F04/F16/F28は今回の依頼で事前指定された視覚的変動の表示例。元contact-frequency値をそのまま用い、統計的有意性・rho順位から選ばない。biomarkerやhitとは呼ばない。他の44 featureもB/Eに残る。Cの色尺度は0〜1でB/Eと異なる。
- **D — External biological anchor: Nishino 2008**: diazepam/triazolamだけの2 mg/kg・5 mg/kg経口投与、15/30/60/90分のrotarod failure fraction。元の16行をそのまま使用し、doseごとに表示。mouse/ICR/male、元表の各dose/time条件n=10 mice。既存fractionは前runでfailure count / sample sizeとして作成済みであり、本作業で再計算・再正規化しない。n=2 drugsのcontext-matched comparisonに限り、external anchorはmediatorではない。相関や因果推論は行わず、dose/time点を独立薬剤と数えない。
- **E — Context phenotypes — not primary wobbling endpoint**: C3 Somnolence、C3 Sedation、C4 Dizziness、C4 Vertigo、C2 Fallの順。Bとは別のheatmapで、同じ48 feature・既存rho・n_effectiveを示す。primary C1/C5と混ぜない。C6は本図に表示しないが元Y snapshotには保持する。

「薬剤間で、あるPLIF featureとclinical balance/motor termが同方向または逆方向に変動する」という対応の記述に限る。残基機序、予測性能、fallの原因、rotarod-mediated clinical effectとは解釈しない。複数C1 termの符号一致は関連した用語間の記述的傾向であり、独立再現や有意性の証拠ではない。FAERS logRORは報告disproportionalityで、incidenceやrisk ratioではない。alpha2 provisional construct、n=10、多重な探索的比較等の既存制約は残る。

## Exact sources and mappings

Source run: `runs/20260920_clinical_bridge_hypothesis_poc/`, introduced at commit `a39e7668ec31eae0e5fed6a0488b97d0815576f4`.

| Input relative to source run | Local byte-identical snapshot | Use |
|---|---|---|
| `raw/lineage/structural_fingerprint.csv` | `data/structural_fingerprint.csv` | Full frozen X, original feature identity/order; C values |
| `tables/clinical_fingerprint_10drug.csv` | `data/clinical_fingerprint_10drug.csv` | Full frozen Y identity/order and NA audit; not recalculated |
| `tables/structural_feature_clinical_correlations.csv` | `data/structural_feature_clinical_correlations.csv` | All 672 archived rho/n values; B/E extraction and descriptive C1 summary |
| `tables/invivo_anchor_contexts.csv` | `data/invivo_anchor_contexts.csv` | D: original `value` by `drug`, `dose`, `observation_time` |
| `STRUCTURAL_RESULTS_FROZEN.json` | same name under `data/` | X SHA256 verification |
| `CLINICAL_RESULTS_FROZEN.json` | same name under `data/` | Y SHA256 verification |
| `config/clinical_phenotype_config.json` | `data/clinical_phenotype_config.json` | Original clinical definitions |
| `config/ANALYSIS_CONTRACT_source.md` | `data/ANALYSIS_CONTRACT_source.md` | Analysis contract provenance |

`input_manifest.json` records all exact source paths and SHA256. Original X originates from `runs/20260919_0724_10drug_standardized_redocking/`. User instructions are saved in `USER_REQUEST.md`. This run is a visualization/summary extension of the same frozen X↔Y comparison, not a revision to scientific feature definitions. Historical runs remain unchanged.

X SHA256: `c986ba328a01fd4bc6909ca24e4724d7eea5971e5a5ea6ae0fc2faa3147b2f10`  
Y SHA256: `8475dd9e413a97075820b93531550ffa8be49dba591834dab3f4e5ccdcb01ea2`

`structural_feature` maps to `feature_label`; its 1-based X column index maps to F01–F48. `phenotype_axis`→`clinical_axis`, `phenotype_term`→`clinical_term`, `rho`→`spearman_rho`; `n_effective` is copied unchanged. `tables/feature_display_key.csv` maps every F identifier to the full source label. Alpha1 is F01–F21; alpha2 is F22–F48.

Drug display order is frozen Y row order: diazepam, alprazolam, triazolam, zolpidem, lorazepam, clonazepam, midazolam, temazepam, zopiclone, zaleplon. See `tables/drug_display_order.csv`. No clustering or result-dependent sorting is used.

## Descriptive tables and missingness

- `tables/wobbling_plif_clinical_correspondence.csv`: 288 rows = all 48 features × 6 C1/C5 terms. Requested six columns: feature_id, feature_label, clinical_axis, clinical_term, spearman_rho, n_effective.
- `tables/context_plif_clinical_correspondence.csv`: separate 240 context cells, same schema, never combined into the primary endpoint.
- `tables/c1_direction_consistency_summary.csv`: 48 rows in F01–F48 order. For the four C1 rho values only, counts of estimable / positive (>0) / negative (<0) / exactly zero, median, minimum and maximum. This is the **only new numerical summarization**, explicitly requested; no individual correlation is recomputed. NA values are excluded from those descriptive summaries, zero remains estimable but neither positive nor negative. If no C1 rho is estimable, counts are zero and median/min/max are NA. No significance score, p-value ranking, multiple-testing hit call or threshold-based feature retention.
- `tables/representative_frozen_plif_values.csv`: unchanged original values for the four user-designated C examples, in the same 10-drug order.

Original Y has six NA cells from its existing low-count rule. They are not reinterpreted or filled. Undefined archived rho remains NA even when n_effective is available (e.g. constant X); n alone is not evidence of an estimable correlation. All 48 features remain visible, including fully undefined rows.

## Reproduce and outputs

Requires Python + numpy + pandas + matplotlib. From repository root:

```bash
python runs/20260920_wobbling_focused_correspondence/scripts/plot_wobbling_focused_correspondence.py
```

The run directory can also be regenerated from bundled data snapshots without network access. The script verifies saved input hashes and frozen X/Y hashes, unique identities, source dimensions, term mappings and anchor contexts. Each exported rho/n cell is checked against its archive. PNG is 200 dpi (4800×5600); SVG heatmaps use vector cells and preserve text. Outputs:

- `figures/wobbling_focused_plif_clinical_correspondence.png`
- `figures/wobbling_focused_plif_clinical_correspondence.svg`
- tables listed above, `input_manifest.json`, `QC_REPORT.json`

[Nishino2008](https://doi.org/10.1254/jphs.08107FP): previously audited Table 1 values; no new extraction in this task.
