# IDを固定したChEMBLデータ取得

Python 3.10以降の標準ライブラリで実行できる。既存解析の追加依存（Vina、RDKit、PLIP）は不要。

```bash
python3 -m src.acquire_chembl
# 別registry等を参照するJSON設定で拡張
python3 -m src.acquire_chembl --config config/chembl_acquisition.json
# 保存済み応答だけを使い、新規runに再処理（ネットワーク不要）
python3 -m src.acquire_chembl --replay-run runs/<run_id>
# raw応答とmanifest記録済み成果物のハッシュ確認（読み取り専用）
python3 -m src.acquire_chembl --verify-run runs/<run_id>
# テスト（macOSでもheadless backendを明示）
MPLBACKEND=Agg python3 -m pytest -q
```

各実行は `runs/YYYYMMDD_HHMM_chembl_acquisition[_N]/` を作り、過去runを上書きしない。`runs/LATEST_RUN.txt` のみ現runへ更新される。異常終了も `FAILED` として保存されるため、ポインタだけではなくmanifestの状態を確認する。replayでは元runの設定・registry snapshotを用い、元の取得日時を保持する。再取得日時と取り違えない。

## 化合物同定と追加手順

1. `config/compound_registry.csv` に既存と重複しない `drug_id` と **明示的な `molecule_chembl_id`** を追加する。名称だけの行や重複IDはエラーとする。
2. 可能なら独立に確認したfull `expected_standard_inchikey`、識別子の出典 `identity_source`、根拠commit/ノートを記録する。新IDの選定は構造・立体・塩を確認して行う。検索結果から名前だけで選択する機能は設けていない。
3. IDに対応するmoleculeレコードを取得する。expected InChIKeyがある場合は完全一致を必須とし、不一致/欠落時は活性取得を停止してQC failureとする。expected keyなしの場合は `explicit_id_only_structure_review` と明記し、構造照合済みとは扱わない。
4. parent/active ID、chirality、full InChIKey、canonical SMILES、複数componentの有無を保持する。親化合物・塩・立体異性体のレコードを自動統合せず、親IDへqueryを拡張しない。
5. QC後に別フェーズでIFP作成候補を判断する。`ifp_available` は入力registryの注記であり、活性取得から推論しない。

## 検索範囲

公式ChEMBL APIのmolecule、target、activity、assay、document、statusを使う。対象targetは設定された名称検索語（GABA-A関連表記、aminobutyric acid receptor、benzodiazepine）と明示target IDの和集合。APIが返した候補と除外理由は `tables/target_scope.csv` に保持する。

関連蛋白・peripheral/translocator候補を除外し、中枢・末梢をまとめた混合targetは要確認区分として保持する。GABA-A、rho、benzodiazepine名称候補を別scopeとして保持する。GABA-A関連というだけではBZD-site活性と確定しない。名称検索範囲外、別の標的表現、親/塩IDの活動は取得漏れになり得る。**ChEMBL全体のBZD-siteデータを網羅したという主張はしない。**

activityはcompound ID × target ID batchで取得し、activity type、単位、種を限定しない。全ページを巡回し、総件数・ページ進行・返却molecule/target IDを検証する。API停止、ページ不足、異常schemaは「0件」にせず失敗する。HTTP 429/一時エラー・タイムアウトには限定回数のretryを行い履歴を保存する。失敗時は新しいrunで再試行する。取得中にChEMBL releaseが変わった場合も失敗とする。

## 保存形式とprovenance

| 場所 | 内容 |
|---|---|
| `config/` | 設定、registry、取得コード/run managerのsnapshot |
| `raw/responses/*.json.gz` | 成功応答のJSON本文を変更せずgzip圧縮。`gzip -dc`で原文を復元可能 |
| `raw/requests.jsonl` | URL（全query付き）、GET、取得UTC時刻、HTTP status/headers、retry、原文SHA-256、圧縮ファイルSHA-256、replay元 |
| `processed/compounds.{csv,jsonl}` | ID、SMILES、InChIKey、hierarchy、構造照合status、取得元 |
| `processed/targets.{csv,jsonl}` | 全targetレコード、component/accession/synonym/xref情報 |
| `processed/assays.{csv,jsonl}` | 全assayレコード、description、parameters、organism、cell/tissue、variant、confidence等 |
| `processed/documents.{csv,jsonl}` | 文献ID、DOI、PubMed ID等（原レコードにある項目） |
| `processed/activities.{csv,jsonl}` | 全activityフィールド＋assay・文献詳細・targetへの参照。raw/standard value・type・unit・relationを別列で保持 |
| `tables/` | 化合物QC、行単位QC、条件別coverage、subunit要確認表、件数summary |
| `run_manifest.json` | 基準Git commit、実行コードhash、input snapshot hash、出力hash、lineage、版、状態 |

processedの各行は `source_request_id`、`source_file`、`source_json_pointer`、`source_sha256`、`retrieved_at_utc` で元応答へ戻れる。joinされたassay/target/documentにも個別のprovenanceがある。targetの大きなxrefはtargets表に一度だけ格納する。JSONLを型・nullを保持する基準形式とし、CSVのdict/list列はJSON文字列、欠損は空欄とする。CSVではnullと原文の空文字の区別が失われるため、必要時はJSONL/原応答を使用する。

`chembl_version` はChEMBLのsource標準化版であり、この処理で構造・単位・活性値を独自に標準化したという意味ではない。manifestのGit commitは実行開始時の基準commit。未commitコードもsnapshotとhashで識別できる。

## QCと解釈

- moleculeのID/構造整合、SMILES・InChIKey欠損、parent差、複数component、chiralityを記録する。SMILESの化学的妥当性をRDKitで検証する処理は今回は含めない。
- activity IDの反復、ChEMBL `potential_duplicate`、`data_validity_comment` を保持する。重複候補を自動削除しない。
- numeric value、unit、relation等の欠損、非数値・非有限、censoringをフラグ化する。`<`/`>`を等号へ変えず、log/pChEMBL換算・平均化・imputationを行わない。unit欠損には無次元指標もあり、すべてをエラーと断定しない。
- assayとactivityのtarget IDを照合し、assay/documentの参照欠損を検査する。target名のsubunitとassay本文のsubunitは **文字列から抽出したreview hint** として分離する。
- 単一subunit targetにβ/γの説明が加わる場合も `text_discrepancy_review` に含まれる。同じfamilyの番号が異なる場合は追加で `subunit_conflict_review`。複数受容体を比較する記述もあり、flagだけで誤注釈と断定・修正しない。原文・原著の確認が必要。
- target component一覧から、実際のassayで使用された組成・stoichiometryを推定しない。構造化assay parameterが空の場合、自由記載descriptionを保持し、pH・温度・濃度などを推測しない。
- `COMPLETE_WITH_REVIEW_FLAGS` は記録された検索範囲の取得完了であり、学習可能性や比較可能性の合格判定ではない。identity/構造欠損、参照欠落、activity ID反復、assay target ID矛盾は `PARTIAL_QC_FAILURE`（終了コード2）。取得自体の失敗は `FAILED`（非ゼロ終了）。

将来のMLではcompound構造、target・subunit、species、assay、document、endpoint・unit・relationを維持して採用条件を決め、同一compound/出典に由来する漏洩を避ける必要がある。今回の全activity表をそのままtraining labelとしない。

## 参照した公式資料

- [ChEMBL Data Web Services](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services)
- [ChEMBL REST API specification](https://www.ebi.ac.uk/chembl/api/data/docs)

実取得のURL・版・応答は各runのraw領域を参照する。
