# Archived Vina docking-score layer

This run recovers archived Vina docking scores and establishes a separate descriptive layer. No docking, PLIF calculation, FAERS retrieval, frozen X/Y modification, feature selection, or numerical PLIF–Vina integration was performed. The task-specific extension to the historical analysis contract is documented in USER_REQUEST.md: seed summaries and receptor-specific docking-score/C1 correspondence are authorized; the historical contract and results remain unchanged.

## Source identity and exact provenance

Only `runs/20260919_0724_10drug_standardized_redocking/raw/archives/{drug}.tar.gz` was read for docking. Each archive contains `{drug}/{receptor}/{seed}/poses.pdbqt`, where receptor is alpha1 or alpha2 and seed is 2026091901–2026091905. Drugs, in frozen order: diazepam, alprazolam, triazolam, zolpidem, lorazepam, clonazepam, midazolam, temazepam, zopiclone, zaleplon.

This source is the standardized 10-drug run underlying frozen PLIF X, rather than the historical 4-drug PoC. `config/docking_run_manifest_source.csv` identifies every drug/receptor/seed output; `config/raw_archive_manifest_source.csv` links its archived member. All 100 PDBQT byte hashes match both manifests. The 10 archive hashes match the standardized run's `run_manifest.json` artifact hashes. `data/archive_manifest.csv` and `data/pdbqt_file_manifest.csv` enumerate exact source paths, expected/observed hashes, pose counts and snapshot paths. Unmodified extracted bytes are in `data/pdbqt/{drug}/{receptor}/{seed}/poses.pdbqt`. A source path using `archive::member` uniquely identifies the original archived file.

Frozen inputs come from `runs/20260920_clinical_bridge_hypothesis_poc/raw/lineage/structural_fingerprint.csv` and `tables/clinical_fingerprint_10drug.csv`. The former contains 21 alpha1 and 27 alpha2 columns. Byte-identical snapshots and freeze metadata are included; `config/frozen_input_manifest.json` records provenance. Display order is copied from `runs/20260920_wobbling_focused_correspondence/tables/drug_display_order.csv`.

## Extraction and summaries

The parser preserves every REMARK VINA RESULT record, using MODEL as pose rank and recording score, both RMSD bounds, line number and raw text. It flags malformed/nonfinite numbers, invalid RMSD bounds, missing/duplicate MODEL records, file/seed duplication or absence, unexpected pose counts and provenance mismatch. Malformed records are retained in the raw table and flagged; unreliable seed outputs are not silently treated as valid. The parser builds on the archived-score recovery script saved in `config/original_recovery_script.py`.

Recovered: **100 files, 883 poses, 100 seed-best rows, 20 drug/receptor summaries**. All conditions have five seeds. File pose counts: 91 × 9, 3 × 8, 4 × 7, 2 × 6; all match their source manifest (up to 9 poses is expected). No missing seeds, malformed records or provenance issues were observed. `data/parsing_issues.csv` and `data/missing_seeds.csv` are header-only, not omitted.

Within each seed, select the minimum Vina docking score; all 100 minima agree with rank 1. Preserve all poses and all five seed minima. Primary drug/receptor value is the median of seed minima, not the global minimum. Mean, sample SD (ddof=1), linear-interpolated IQR (Q75−Q25), minimum, maximum and n_seeds_available are secondary summaries. `vina_favorability = -median_best_seed_score` is an explicit sign reversal; raw scores remain present. Docking tables are hashed in VINA_LAYER_FROZEN.json before clinical Y is loaded.

## Descriptive correspondence

Compare drugs with Spearman rho, separately for alpha1 and alpha2. C1 terms are fixed: Ataxia, Balance disorder, Coordination abnormal, Gait disturbance. Use pairwise complete finite values without imputation; require n >= 5 and nonconstant inputs, otherwise return NA. n_effective is 9, 10, 8, 10 respectively. Frozen Y is logROR; its existing low-count masks are unchanged (zaleplon Ataxia; triazolam and zaleplon Coordination abnormal). Raw-score and favorability correlations contain identical rank information and opposite signs. There are no p-values, cutoffs, significance scores, or hit calls.

| Receptor | Ataxia (n=9) | Balance disorder (n=10) | Coordination abnormal (n=8) | Gait disturbance (n=10) |
|---|---:|---:|---:|---:|
| alpha1 favorability rho | -0.366667 | -0.660606 | -0.357143 | -0.709091 |
| alpha2 favorability rho | -0.200000 | +0.260606 | 0.000000 | -0.006061 |

Alpha1 versus alpha2 raw median-score rho = -0.151515 (n=10). This describes ordering across the ten drugs and does not establish receptor equivalence. Seed robustness retains every seed, including the higher zopiclone/alpha2 score around -7.33 while the other four are around -8.1; no outlier exclusion is applied.

## Outputs

- `data/vina_scores_all_poses.csv`: all 883 poses, original source references and QC fields.
- `data/vina_scores_seed_best.csv`: 100 independent seed minima, original and requested best-score column names.
- `data/vina_scores_drug_receptor_summary.csv`: 20 rows, all five seed scores, median and other summaries, derived favorability.
- `data/vina_alpha1_alpha2_correspondence.csv`: descriptive cross-block rho.
- `data/vina_clinical_c1_correspondence.csv`: eight receptor/term rows with raw and favorability rho, n, included/excluded drug IDs and status.
- `data/vina_clinical_c1_drug_values.csv`: 80 input pairs and inclusion masks for audit.
- `data/structural_multiview_index.csv`: ten drug rows with two PLIF vector identifiers, separate median scores and favorabilities. Identifiers resolve through `data/plif_vector_reference_manifest.json` to exact frozen-file hashes, row keys and ordered feature columns. This index joins references, not numeric features.
- `figures/vina_score_alpha1_alpha2_heatmap.png/.svg`: frozen drug order, raw seed-median score and n_seeds.
- `figures/vina_score_seed_robustness.png/.svg`: all seeds in separate receptor panels; black bars are medians. Panel y-scales differ.
- `figures/vina_score_alpha1_vs_alpha2.png/.svg`: labeled drug scatter with descriptive rho.
- `figures/vina_clinical_c1_correspondence.png/.svg`: favorability rho colors and n_effective cell annotations.
- `QC_REPORT.json`, provenance manifests, configuration snapshots, `VINA_LAYER_FROZEN.json`, and scripts support reproducibility.

## Interpretation limits

Vina docking score is a computational docking-score proxy, not experimentally measured binding affinity or true binding free energy. More negative scores are more favorable within this scoring setup. Seeds are computational search repeats, not biological replicates. Alpha1 and alpha2 remain separate; alpha2 retains the provisional construct and beta3-related structural limitations documented in the source contract.

Observational unit = drug. Clinical Y retains its frozen all-role FAERS reporting limitations, missing cells, lack of patient-level covariate adjustment and reporting confounding. These ten-drug correlations are descriptive QC/hypothesis generation only, without causal, mediation, mechanistic or predictive claims. Association/correspondence only; not causation or mediation. No PLIF weighting, composite score, clinical optimization, residue selection or PLIF reranking is performed.

## Reproduction and validation

From repository root, with Python plus numpy, pandas, scipy and matplotlib:

```sh
python runs/20260920_vina_score_layer/scripts/recover_vina_score_layer.py
python runs/20260920_vina_score_layer/scripts/test_recovery.py
```

The script requires the standardized archived source run and frozen-input source runs at the paths above. It checks existing snapshots before reuse. Re-execution regenerates only this run's derived files and refreshes the freeze timestamp. `config/software_versions.json` records the environment used.

Ten tests pass: six parser edge cases, independent comparison of all archived numerical records and seed summaries, sample-SD/favorability arithmetic, independent rank-Pearson verification of all C1 Spearman values, and frozen docking hashes plus all twenty PLIF vector references. All four PNG figures were visually inspected for labels, clipping and coverage; SVG counterparts were also generated. Historical tracked files remain unchanged.
