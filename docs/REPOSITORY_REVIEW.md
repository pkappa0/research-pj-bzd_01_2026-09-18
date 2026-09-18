# 既存解析のレビュー（2026-09-18）

対象は `pkappa0/research-pj-bzd_01_2026-09-17`、確認時点のmainは `dc2c07ded7a3655ff3176ad1db40007b1a90fa83`。Git履歴はこの1 commitのみで、過去の解析工程はcommit履歴ではなく `runs/` のmanifest・FINAL_REPORT・QC_REPORTから追跡した。既存の全追跡CSVの行数・列数・SHA-256は `baseline_csv_inventory.csv` に記録した。

## 現状

- **構造IFPは4剤**（diazepam、alprazolam、triazolam、zolpidem）。α1は6HUP、α2は9CTJ由来の暫定的な局所construct。4剤 × 2系 × 9 pose = 72 docking行、PLIPは339相互作用行。
- `20260917_1604_multireceptor_feature` では177行がmapping未解決だった。後続 `20260917_1629_residue_hypothesis` でγ2鎖を補足し、`20260917_1706_structure_mouse_bridge` で受容体を含む鎖・残基キーを使用した結果、339行すべての対応が解決された。初期reportの「未解決177行」を最新状態として扱わない。
- 最新の既存run `20260917_1926_residue_class_heatmap` は25行の候補残基・interaction特徴を再表示したもの。raw frequencyと行内中心化値を分離し、後者は図示用途のみ。
- phenotype bridgeはマウスの文献由来44行、数値primary matrixは8行。同一条件の4剤数値比較は成立せず、alprazolam・zolpidem等の欠損を補完していない。これらは本フェーズでは再抽出・再検証していない。
- α2暫定construct、box transfer、alprazolamのβ2薬理対β3 docking背景は解釈上の制約。mapping解決だけで実験条件が揃うわけではない。ML学習可能性・予測精度は未検証。

## 化合物情報

`data/raw/drug_list.csv` および `data/raw/pharmacology/chembl/molecules.csv` には、既に10剤のID・SMILES・InChIKeyがある。薬理raw CSVは434行、39 target。今回の取得対象はこの既存10剤を明示IDで固定し、IFP未作成の6剤も同じ取得・QC契約で扱えるようにするもの。新しい化合物を名前だけで追加してはいない。

| drug_id | ChEMBL ID | 既存IFP | 既存活性行 |
|---|---|---|---:|
| diazepam | CHEMBL12 | あり | 263 |
| lorazepam | CHEMBL580 | なし | 2 |
| clonazepam | CHEMBL452 | なし | 16 |
| alprazolam | CHEMBL661 | あり | 16 |
| midazolam | CHEMBL655 | なし | 3 |
| temazepam | CHEMBL967 | なし | 3 |
| triazolam | CHEMBL646 | あり | 12 |
| zolpidem | CHEMBL911 | あり | 109 |
| zopiclone | CHEMBL135400 | なし | 4 |
| zaleplon | CHEMBL1521 | なし | 6 |

ID assignmentの出典をregistryに固定し、取得時にはfull InChIKeyも照合する。既存の名前一致による割当履歴を無条件に再採用するのではなく、IDと保存済み構造識別子の整合性を確認する。名称は表示用でありjoinキーではない。既存IFP側のdrug_idはregistryを経由してChEMBL IDへ対応づける。

## 既存コードと出力の注意

- `src/fetch_chembl.py` はpreferred name/synonymの一意exact matchを採用する。元APIページの保存はあるが、同じ出力先への上書き、原ページごとのURL・ハッシュ不足、activity加工表でのassay詳細欠落がある。旧実行コマンドを今回は再実行しない。
- `src/strict_poc.py`、`src/four_drug_poc.py` は構造・薬理の旧工程。`src/multireceptor_poc.py`、`src/residue_hypothesis.py`、`src/structure_mouse_bridge.py` とheatmapスクリプトはrun単位の後続工程。今回はこれらを変更・再実行しない。
- `src/run_manager.py` の新規run生成、lineage、code snapshot、FINAL/QC reportの方針を新パイプラインでも使用する。
- 旧 `outputs/FINAL_REPORT.md` はstrict 3剤段階の報告であり、最新4剤・mapping補正済み状態を表していない。
- 旧READMEの「drug_listがない」は過去時点の記述。現リポジトリではリストも10剤の取得結果も存在する。
- `20260917_1602_base_pipeline` は入力欠損run。新しい日付や `LATEST_RUN.txt` だけで有効な科学的成果物を選ばない。
- `src/fetch_gtopdb.py` は補完機構を持つが、既存reportではAPI key/fallback CSVがなく取得0行。今回はChEMBLの取得を実装・検証し、GtoPdbを取得済みと扱わない。
- 原PDB/PLIP実行環境など、Gitに含まれないものや旧絶対パスを参照する箇所がある。今回のcloneだけで過去の全dockingが再現可能という意味ではない。

## 今回の方針

独立した `src/acquire_chembl.py` と `config/compound_registry.csv` を追加する。旧434行や既存fingerprint、phenotype、過去runは更新しない。新規runで取得完了性と科学的な比較可能性を分け、確認が必要な行を残す。詳細は `DATA_ACQUISITION.md` を参照。
