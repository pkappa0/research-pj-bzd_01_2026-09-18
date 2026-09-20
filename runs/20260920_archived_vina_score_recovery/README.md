# Archived Vina docking score recovery

既存standardized 10-drug runのPDBQTに保存された **Vina docking score** の回収。**ドッキング、PLIF、clinical RORは再実行・再計算していない。** alpha1/alpha2を別条件として保持し、臨床phenotypeによる調整、最適化、rescaling、受容体間の合成scoreを作らない。Vina docking scoreは計算値であり、実験で測定した結合親和性ではない。

## Recovery outcome

- Source run: `runs/20260919_0724_10drug_standardized_redocking/`.
- Parsed 100 PDBQT outputs = 10 drugs × 2 receptor blocks × 5 computational search seeds.
- Preserved **883 pose records**, including all non-best poses.
- Output pose counts: 91 files have 9 poses, 3 have 8, 4 have 7, 2 have 6. Counts match the original docking manifest. Nine is an upper bound, not a reason to fill in absent poses.
- All 100 seed-level rank-1 scores equal the minimum score in their respective file.
- All 20 drug/receptor conditions retain 5/5 seed-best scores.
- Missing seeds: 0. Malformed result records: 0. Other parsing/provenance flags: 0.
- All 10 archive hashes and all 100 PDBQT hashes match archived provenance. Files are validated against both the docking output manifest and the archive member manifest.

## Exact inputs / provenance

Docking outputs are tar members:

`runs/20260919_0724_10drug_standardized_redocking/raw/archives/{drug}.tar.gz::{drug}/{alpha1|alpha2}/{seed}/poses.pdbqt`

Every available PDBQT member in these ten archives is scanned. Receptor/ligand preparation PDBQT inputs under `raw/receptors` and `raw/ligands` are not docking outputs and are not included. The 100 exact parsed paths, archive-member byte SHA256, manifest checks, pose counts and copied snapshot paths are in **`data/pdbqt_file_manifest.csv`**. The ten archive hashes are in `data/archive_manifest.csv`. Each original output is copied byte-for-byte under `data/pdbqt/{drug}/{receptor}/{seed}/poses.pdbqt`.

Original `config/docking_config.json`, `tables/docking_run_manifest.csv`, and `tables/raw_archive_manifest.csv` are preserved as source snapshots under `config/`. `config/input_manifest.json` records paths and hashes. `config/ANALYSIS_CONTRACT_source.md` records the repository analysis contract: compare drugs within fixed receptor blocks; alpha1/alpha2 are not replicates. Original outputs and historical runs are not changed.

Frozen display order is read only from `runs/20260920_wobbling_focused_correspondence/tables/drug_display_order.csv` and snapshotted under `config/`. No clinical values are read. Order: diazepam, alprazolam, triazolam, zolpidem, lorazepam, clonazepam, midazolam, temazepam, zopiclone, zaleplon.

## Parsing and aggregation rules

- Scan every recognizable `REMARK VINA RESULT` line, require `REMARK VINA RESULT:` followed by exactly three finite numeric fields. They are copied as `vina_score_kcal_mol`, `rmsd_lb`, and `rmsd_ub` (RMSD bounds in angstroms). Source `MODEL` integer defines `pose_rank`; the score is not used to renumber poses.
- Long table preserves **every** result record, including malformed recognizable lines if encountered. Additional fields `line_number`, `record_index`, `record_status`, and `raw_record` provide auditability. Unparseable numeric records have NA numeric fields, not zero. Parsed but invalid numeric values are retained with flags.
- Validate model/rank structure, one result per model, RMSD nonnegativity/order, source hashes and original manifest pose counts. Missing MODEL rank, duplicate rank, missing/multiple records, malformed fields, unexpected identity, hash discrepancy and multiple files per seed are explicitly flagged.
- Within each drug/receptor/seed, the best score is the minimum (most negative). Lowest original pose rank breaks exact ties. Rank 1 is independently checked and its value/agreement retained; all rank-1/minimum comparisons agree in this dataset. No other poses are discarded from the raw table.
- Files with structural parsing/provenance problems are retained in the raw table but not summarized. Expected missing seeds have explicit NA rows in the seed table and appear in `data/missing_seeds.csv`. Missing data are never imputed. `data/parsing_issues.csv` retains detailed flags; both flag CSVs contain headers and no issue rows for this run.
- Summary statistics use the five available independent computational search-seed best scores **within each drug and receptor**. Median, mean, minimum, maximum and IQR are in kcal/mol. IQR = Q75 − Q25 using linear interpolation (NumPy/pandas type 7); for five sorted scores this is fourth minus second. No standardization or cross-receptor aggregation is performed.
- These seeds are computational search repeats, not biological replicates. Variability describes this fixed search protocol, not experimental affinity uncertainty. Structural construct/preparation differences (including provisional alpha2) prevent interpreting receptor-block score differences as measured selectivity.

## Outputs

- `data/vina_scores_all_poses.csv`: 883 rows; requested drug, receptor_block, seed, pose_rank, vina_score_kcal_mol, rmsd_lb, rmsd_ub, source_file plus audit fields.
- `data/vina_scores_seed_best.csv`: 100 rows, all five seed-level scores for each drug/receptor; rank-1 agreement and ties included.
- `data/vina_scores_drug_receptor_summary.csv`: 20 rows, requested summary statistics and n_seeds_available, plus n_seeds_expected and a separate score column for each of the five seeds.
- `figures/vina_score_median_qc.png` / `.svg`: frozen drug order × alpha1/alpha2; median and n_seeds in every cell. A shared color key is only a display legend, not score normalization or receptor aggregation.
- `figures/vina_seed_best_robustness.png` / `.svg`: all five best scores per condition, separate receptor panels. Small fixed vertical offsets separate seeds. Horizontal ranges show min/max. No clinical optimization.
- `QC_REPORT.json`, `VALIDATION.txt`, exact-file manifests and issue CSVs.

## Reproduce / verify

Requires Python with numpy, pandas and matplotlib; no Vina executable or API access is needed. Run from repository root:

```bash
python runs/20260920_archived_vina_score_recovery/scripts/recover_vina_scores.py
python runs/20260920_archived_vina_score_recovery/scripts/test_recovery.py
```

The recovery script reads the existing repository archives; standalone copies of parsed PDBQTs are additionally provided for audit. Tests cover valid/non-best preservation, malformed numeric rows, missing MODEL/result records, invalid RMSD, duplicate result records, and independent validation of every archived RESULT field and all summary arithmetic. Tests do not run docking.
