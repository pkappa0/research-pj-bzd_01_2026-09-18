# GABA-A Interaction Fingerprint PoC

GABA-A receptor docking/PLIP結果を入力として、残基レベルのInteraction Fingerprintを作り、薬理活性・in vivo/clinical phenotypeを同じ薬剤順で並べて探索的に確認するための新規PoCリポジトリです。ドッキング計算そのものは実行せず、既存のVina/PLIP出力を後から配置して解析します。実データがない状態でも処理を止めず、不足入力を `results/tables/qc_summary.csv` に出力します。

## 実行

```bash
python -m pip install -r requirements.txt
python run_pipeline.py --config config/poc.yaml
pytest -q
```

`tests/fixtures/` のデータは動作確認専用のsynthetic dataです。研究結果・薬理学的主張・実データの代用として扱わないでください。通常の解析対象は `data/raw/` に置きます。

## ディレクトリ

```text
data/raw/interactions/interactions.csv       # PLIP等のcontact表
data/raw/residue_mapping.csv                # subtype間の対応
data/raw/pharmacology/pharmacology.csv
data/raw/phenotype/phenotype.csv
data/processed/                              # 自動生成CSV
results/figures/                             # 自動生成PNG
results/tables/                              # QC・探索統計
src/                                         # 処理モジュール
tests/fixtures/                              # synthetic専用
```

## 実データ入力テンプレート

テンプレートは `data/raw/**/` にあります。ファイル名をテンプレートから実ファイル名へ変更するか、`config/poc.yaml` のパスを変更してください。ID・名称・指標の欠損を0で埋めないでください。

### 1. `data/raw/interactions/interactions.csv`

必須列は `drug_id` (string)、`receptor` (string: 例 `alpha1`, `alpha2`)、`pose_id` (string)、`original_residue` (string: 例 `TYR160`)、`interaction_type` (string) です。`interaction_type` はPLIP等の実出力にある名称をそのまま記録してください（例 `hydrogen_bond`, `hydrophobic`, `pi_stacking`, `pi_cation`, `salt_bridge`, `halogen_bond`）。

任意列は `drug_name_raw`、`pdb_id`、`chain`、`residue_name`、`residue_number`、`source_file`、`score`、`confidence` です。1行は1つのdrug × receptor × pose × residue × interaction typeです。同じpose内の重複行は読み込み時に除去します。

### 2. `data/raw/residue_mapping.csv`

必須列は `receptor`、`original_residue`、`common_position` (string) です。`common_position` はPDB番号ではなく、subtype間で検証済みの共通位置ID（例 `POS_023`）を用意してください。任意列は `mapping_confidence`、`mapping_note`、`alignment_source` です。対応が曖昧な残基は高い確信度として登録せず、noteに理由を書いてください。mappingがないinteractionは `UNMAPPED_<receptor>_<residue>` として保持され、QCに件数が出ます。

### 3. `data/raw/pharmacology/pharmacology.csv`

必須列は `drug_id` (string) です。`drug_name_raw` と、測定指標ごとの列（例 `ki_alpha1_nM`, `kd_alpha1_nM`, `ic50_alpha1_nM`, `ec50_alpha1_nM`, `selectivity_alpha1_alpha2`）を任意列として追加してください。Ki/Kd/IC50/EC50を同じ列や同じ尺度に混ぜないでください。単位を列名または別の明示列に記載し、測定法・文献識別子などの出典列も推奨します。値はraw値のまま保存し、log変換が必要な場合は別列を追加してください。

### 4. `data/raw/phenotype/phenotype.csv`

必須列は `drug_id` です。入力に存在する表現型だけを列として追加し、連続値・binary・categoricalを区別できるよう `phenotype_scale` 等の任意列を付けてください。例は `sedation`, `ataxia`, `amnesia`, `anxiolytic_effect`, `muscle_relaxation` ですが、実際にデータがある項目だけを使用します。文献値を推定して補完しないでください。

## 出力

`data/processed/` は `residue_mapping.csv`、`fingerprint_binary.csv`、`fingerprint_frequency.csv`、`fingerprint_similarity.csv`、`pharmacology_processed.csv`、`phenotype_processed.csv`、`integrated_dataset.csv` を生成します。binaryは `common_position × interaction_type` の存在、frequencyは全pose中の出現割合です。binary similarityはJaccard、frequencyの距離はEuclideanです。Jaccard距離にaverage linkageを適用し、そのdrug orderを全heatmapに固定します。

`results/figures/` にはbinary/frequency/similarity/pharmacology/phenotype/integrated heatmapを生成します。`results/tables/qc_summary.csv` は各入力の `loaded` / `missing` / `invalid`、行数、pose数、mapping未解決行を明示します。`exploratory_correlations.csv` は十分なdrug数と完全な数値pairがある場合だけSpearman探索相関を計算し、因果関係は解釈しません。

## 今回あなたが用意する実データ

1. `data/raw/interactions/interactions.csv` — 全drug、受容体、poseについてのPLIP/contact表。
2. `data/raw/residue_mapping.csv` — alpha1/alpha2の残基対応と、曖昧性・alignment出典。
3. `data/raw/pharmacology/pharmacology.csv` — drugごとのKi/Kd/IC50/EC50/selectivity（指標・単位を分離）。
4. `data/raw/phenotype/phenotype.csv` — 入手済みのin vivo/clinical phenotypeと尺度・出典。
5. 可能なら元のPLIP report、PDB/mmCIF、Vina score、alignmentファイル — CSVの出典追跡とQCに使用。

これらを配置後、`python run_pipeline.py --config config/poc.yaml` を実行してください。入力が不足していても生成物とQCは出力され、存在する範囲で処理が進みます。

## ChEMBL pharmacology取得

対象薬剤は `data/raw/drug_list.csv` に明示してください。必須列は `drug_id` と `input_drug_name` です。テンプレートは `data/raw/drug_list.template.csv` です。synthetic fixtureやREADME中の例は対象リストとして自動使用しません。

```bash
python fetch_chembl.py --config config/poc.yaml
```

取得スクリプトはChEMBL Web Servicesの分子検索・分子詳細・target検索・activity検索を使います。preferred nameまたはsynonymの一意なexact matchだけを自動採用し、fuzzy候補・複数候補は `results/tables/chembl_qc.csv` とraw JSONに残して保留します。塩/多成分SMILES、molecule hierarchy、chirality、target type、subunit-specific/compositeの分類を保持します。Ki/Kd/IC50/EC50はactivity行ごとに保持し、単位変換・統合はしません。

取得物は次の通りです。

- `data/raw/pharmacology/chembl/chembl_api_metadata.json`: 取得UTC時刻、API状態、ChEMBL version/release
- `data/raw/pharmacology/chembl/api_pages/`: APIレスポンスページの原JSON
- `data/raw/pharmacology/chembl/molecule_raw.jsonl`, `activity_raw.jsonl`: 入力drug IDを付加した原レコード
- `data/raw/pharmacology/chembl/molecules.csv`, `activities.csv`: raw抽出表
- `data/raw/pharmacology/chembl/target_candidates.csv`: GABA-A関連target候補。subunit-specificとcompositeを別列で保持
- `data/processed/pharmacology.csv`: 解析用processed activity表
- `results/tables/chembl_qc.csv`: missing、未解決名、候補数、activity件数

現在の作業フォルダには `data/raw/drug_list.csv` がないため、今回の実行ではAPI接続とtarget候補・ChEMBL versionの取得だけを行い、薬剤activity行は生成していません。薬剤リストを配置した後に同じコマンドを再実行してください。

## Strict 3-drug mechanistic PoC

既存のChEMBL raw activityを変更せず、`diazepam`、`triazolam`、`zolpidem`だけを対象に、α1β3γ2（CHEMBL2094121）とα2β3γ2（CHEMBL2094130）のstrict解析を実行する独立スクリプトがあります。

```bash
python src/strict_poc.py
```

この実行は、同一document・Ki・unitを優先したpair-level pharmacology、ChEMBL canonical SMILESからのRDKit 3D conformer、既存RCSB PDBを使ったVina、vendored PLIP、sequence alignmentを経由した残基mapping、binary/frequency IFP、heatmapと統合figureを生成します。6HUPは直接のhuman α1β3γ2構造です。exact human α2β3γ2 full pentamerにBZD共結晶ligandを持つ候補を確認できなかったため、α2側はnative mixed 9CTJからβ3/α2/γ2局所サイトを抽出した暫定constructです。`outputs/receptor_selection_report.md` と `outputs/QC_REPORT.md` に構造・box・mapping・PDBQTの制約を明記し、α2の結果をclean comparisonとは扱いません。

strict成果物は `outputs/` に出力されます。

- `strict_pharmacology_summary.csv`: raw Ki pairごとのα2/α1 ratio。平均値や代表値は作りません。
- `docking_results.csv`: drug × receptor × poseのVina score、seed、条件、PDB/PDBQTパス。
- `residue_mapping.csv`: global sequence alignment + PDB residue lookupのmapping監査表。
- `fingerprint_binary.csv`, `fingerprint_frequency.csv`: pose単位の相互作用fingerprint。
- `alpha1_alpha2_similarity.csv`: drugごとのJaccardと変化feature数。
- `figures/fingerprint_binary_heatmap.png`, `figures/fingerprint_frequency_heatmap.png`, `figures/ifp_vs_pharmacology_integrated.png`。
- `receptor_selection_report.md`, `QC_REPORT.md`, `FINAL_REPORT.md`。

Vinaがない環境では`--skip-docking`で準備/QCだけを確認できます。PLIPは`work/plip_vendor/`に隔離した実行環境を使います。`data/raw/strict3/`以下のPDB、ligand、receptor PDBQT、Vina pose、PLIP reportは追跡用raw出力として保持します。

## Four-drug Tier-1 extension

既存strict 3-drug成果物を保持したまま、ChEMBLの既存tier分類でTier 1となる `diazepam`、`alprazolam`、`triazolam`、`zolpidem`を対象にアルプラゾラムを追加実行します。

```bash
python four_drug_poc.py
```

このスクリプトは、Tier 1 raw pairと同一metric/unit/reference/assay-context内のmedian summary、既存receptor/box/Vina条件でのalprazolam docking（9 pose × 2 receptor）、既存PLIP、残基mapping、4-drug binary/frequency fingerprint、similarity、differential feature表、heatmapおよび統合figureを生成します。既存の `outputs/docking_results.csv`、`outputs/plip_interactions.csv`、`outputs/fingerprint_*.csv`、`outputs/alpha1_alpha2_similarity.csv` は上書きせず、4剤版は `results/tables/*_4drug.csv` と `data/processed/` に分離して保存します。

主な4剤版出力は `data/processed/primary_tier1_pharmacology_raw.csv`、`data/processed/primary_tier1_pharmacology_summary.csv`、`data/processed/fingerprint_binary.csv`、`data/processed/fingerprint_frequency.csv`、`results/tables/primary_tier1_qc.csv`、`results/tables/alpha1_alpha2_similarity_4drug.csv`、`results/tables/differential_features_by_drug.csv`、`results/tables/drug_specific_features.csv`、`results/tables/ifp_vs_pharmacology_4drug.csv`、`results/figures/*_4drug.png`、`outputs/FOUR_DRUG_POC_REPORT.md`です。最終的な解釈境界、未解決mapping、alpha2 provisional construct、box transfer、alprazolamのβ2薬理対β3 docking組成差はレポートとQC表に記録しています。GtoPdbはこの4剤primaryには統合していません。

## Multi-receptor fingerprint feature-engineering PoC

### Immutable run-directory policy

All newly launched phases use one immutable directory per execution under `runs/<run_id>/`. The run ID is generated from the local timezone as `YYYYMMDD_HHMM_<short_phase_name>`; a numeric suffix is added if the same minute and phase already exist. The directory contains the requested `config/`, `data/`, `raw/`, `processed/`, `results/`, `tables/`, `figures/`, `logs/`, and `reports/` areas, plus root-level `FINAL_REPORT.md`, `QC_REPORT.md`, and `run_manifest.json`. Empty `data/raw`, `data/processed`, `results/tables`, and `results/figures` compatibility directories are also created for downstream phases.

The manifest records timestamps, phase/purpose, parent and source run lineage, drugs/receptors, docking parameters, copied input files, generated output files, a git commit when the workspace is versioned, a code fingerprint, notes, and QC status. `runs/LATEST_RUN.txt` contains only the most recent run ID. A new run is always created for a correction, parameter change, or re-analysis; previous run directories are treated as read-only and are never overwritten.

The general pipeline wrapper now enforces this policy:

```bash
python run_pipeline.py --config config/poc.yaml --phase pipeline
```

Missing inputs are copied only when present and are reported in the run's `tables/qc_summary.csv` and `QC_REPORT.md`; the wrapper does not synthesize data. The multireceptor feature phase is likewise versioned by `python multireceptor_poc.py`. The older `outputs/`, `data/processed/`, and `results/` artifacts remain untouched as legacy records; new work should use the run wrappers.

解析の主軸を「同一受容体内の薬剤間比較」に変更する場合は、既存4-drug結果を再利用して次を実行します。

```bash
python multireceptor_poc.py
```

alpha1（6HUP）とalpha2（9CTJ由来provisional local construct）を別々のfeature blockとして、resolved residue mappingだけから frequency、raw interaction count、PLIP distance statisticsを作成します。既存のstrict/4-drug成果物とalpha1-vs-alpha2 secondary解析は上書きしません。出力は `data/processed/alpha1_drug_feature_matrix.csv`、`alpha2_drug_feature_matrix.csv`、`combined_alpha1_alpha2_feature_matrix.csv`、`drug_pharmacology_labels.csv`、`results/tables/alpha1_interdrug_similarity.csv`、`alpha2_interdrug_similarity.csv`、variable-feature表、`unresolved_feature_qc.csv`、`ml_dataset_readiness_qc.csv`、5種類のfeature heatmap、`outputs/MULTIRECEPTOR_FINGERPRINT_REPORT.md`です。

frequency=0はresolved mapping済みで全poseに相互作用がない場合だけ使用し、unresolved mappingはNAとしてprimary matrixから除外します。距離はPLIP raw schemaの `dist`、`dist-d-a`、`centdist`だけを使用し、距離欠損・相互作用欠如を0距離へ変換しません。n=4のため、MLは実行せず、readiness statusを `NOT READY FOR MODEL TRAINING` と出力します。

## Residue-level BZD interface hypothesis phase

残基レベルの仮説検証は、既存multi-receptor runを変更せず、次のコマンドで新しいrun directoryへ出力します。

```bash
python -m src.residue_hypothesis
```

このphaseでは、既存のalpha1 chain D–alpha2 chain D mappingを保持し、6HUPのgamma2 chain Cと9CTJのgamma2 chain Eを一意なglobal sequence alignmentで追加します。TYR58/PHE77を含むchain mapping不足のPLIP行を回収し、`hydrophobic_interaction`と`pi_stack`を統合せずに残します。pi stackingではcentdist、angle、offset、PLIP type P/T、hydrophobic contactではdistを薬剤・受容体・残基別に集計します。主解析は各受容体内の4薬剤比較であり、alpha1-vs-alpha2 global similarityとMLは実行しません。

出力はrun内の`tables/`に保存されます。`recovered_interactions_qc.csv`はraw attributesとsource row IDを保持し、`no_interaction_observed`と`unresolved_after_correction`を区別します。

## GtoPdb pharmacology supplementation

GtoPdb/IUPHARの補完フェーズはstrict 3-drug PoCとは独立して実行します。

```bash
python fetch_gtopdb.py --config config/poc.yaml
```

REST APIを使う場合だけ、API keyをプロセス環境に設定してください。keyはコード・YAML・CSV・ログ・レポートへ保存しません。

```bash
export GTP_API_KEY='your-key'
python fetch_gtopdb.py --config config/poc.yaml
```

`GTP_API_KEY`がない、またはAPIが401/403を返す場合も処理は停止せず、`outputs/GTOPDB_FETCH_REPORT.md` と `results/tables/gtopdb_ligand_match_qc.csv` に理由と再開方法を記録します。APIを使わない場合は、GtoPdb公式downloadのinteraction CSVを次に置いて同じコマンドを再実行できます。

```text
data/raw/gtopdb/gtopdb_interactions.csv
```

fallback CSVでは少なくとも `drug_id`、`input_drug_name`、`target_name` を用意してください。推奨列は `gtopdb_ligand_id`、`gtopdb_ligand_name`、`gtopdb_target_id`、`target_type`、`species`、`interaction_type`、`action`、`affinity_parameter`、`affinity_value`、`affinity_low`、`affinity_high`、`affinity_units`、`pKi`、`pKd`、`pIC50`、`pEC50`、`alpha_subunit`、`beta_subunit`、`gamma_subunit`、`other_subunits`、`reference_id`、`PMID`、`DOI`、`reference_citation`、`assay_id`です。列がない値は補完せず空欄のまま保持します。

GtoPdb ligandはdrug nameだけで確定せず、既存ChEMBLのcanonical SMILESまたはInChIKeyと照合できた候補だけを自動採用します。nameだけ一致する候補、複数候補、chemical identifierがない候補はQCへ出します。

生成物は `data/raw/gtopdb/api/` の原JSON、`data/raw/gtopdb/gtopdb_interactions_raw.csv`、`data/processed/gtopdb_interactions_raw.csv`、`data/processed/gtopdb_interactions_normalized.csv`、`data/processed/pharmacology_master_raw.csv`、`results/tables/gtopdb_ligand_match_qc.csv`、`gtopdb_drug_target_counts.csv`、`gtopdb_alpha1_alpha2_presence.csv`、`gtopdb_receptor_composition.csv`、`gtopdb_comparable_metrics.csv`、`pharmacology_comparable_metrics.csv`、`pharmacology_duplicate_qc.csv`、`pharmacology_tier_summary.csv`、`outputs/GTOPDB_FETCH_REPORT.md`、`outputs/PHARMACOLOGY_INTEGRATION_REPORT.md`です。ChEMBLとGtoPdbのduplicate候補は削除せず、PMID/DOIを第一優先にフラグ付けします。
