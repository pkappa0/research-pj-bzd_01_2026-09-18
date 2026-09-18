# データ取得・QC結果（2026-09-18）

正式run: [`20260918_1117_chembl_acquisition`](../runs/20260918_1117_chembl_acquisition/FINAL_REPORT.md)。ChEMBL_37（APIが報告するrelease date: 2026-05-01）。状態は `COMPLETE_WITH_REVIEW_FLAGS`。取得範囲の完了を示すが、科学的比較・ML用labelの妥当性を保証するものではない。

## 取得結果

- 10剤すべての明示ChEMBL IDとfull InChIKeyが一致。SMILES・InChIKey欠損なし。構造そのもののRDKit検証は未実施。
- target検索scopeは93件、実際に活性を持つtargetは39件。
- 活性434行、assay335件、文献116件、保存されたAPI応答94件。
- 既存434行とactivity ID集合は完全一致。化合物/target/assay/document ID、standard type/relation/unit、standard valueとpChEMBL値も照合して一致した。
- 4剤IFPに加え、lorazepam、clonazepam、midazolam、temazepam、zopiclone、zaleplonの6剤を明示ID・構造照合・assay詳細付きで取得できる。今回は新しいIFPを計算していない。
- target organism内訳はrat 215、human 193、bovine 21、mouse 5行。異なるspeciesやendpointを統合していない。

## QC所見

| 項目 | 活性行数 | 扱い |
|---|---:|---|
| standard value欠損 | 6 | nullのまま保持 |
| standard unit欠損 | 41 | 無次元指標を含み得るため自動除外しない |
| standard relation欠損 | 5 | 推測しない |
| 非等号relation | 45 | 不等号を保持し等号へ置換しない |
| ChEMBL potential_duplicate | 97 | 元flagを保持、自動削除しない |
| ChEMBL data_validity_commentあり | 5 | コメント原文を保持 |
| target subunit記載なし | 146 | 組成を補完しない |
| specific PROTEIN COMPLEXではない | 216 | generic/single targetを別contextとして保持 |
| structured assay parametersなし | 434 | assay descriptionを保持。pH・温度等を生成しない |
| target/assayのsubunit文字列が異なる | 18 | 詳細化・矛盾候補をreview表に保持 |
| 上記のうち同familyの番号差あり | 5 | 原著確認が必要。自動修正しない |

繰返しactivity ID、assay/document参照の欠落、assayとactivityのtarget ID不一致は0件。QC flagは重複して付くので件数を合算しない。

特にactivity `604850`（diazepam）と `595754`（zolpidem）はtarget `CHEMBL2094121` がα1β3γ2である一方、assay本文にはα4β3γ2と記載される。これは既存のprimary α1薬理候補に関係するため、旧labelをそのままMLに用いる前に原著照合が必要。現在の取得処理ではtargetとassayの両方を保存し、どちらが正しいかを推測しない。

さらに `1876214` はtargetのα1とassayのα5、`1701379`・`1701378` はtargetのβ1とassayのβ2に差がある。単一subunit targetにβ/γの説明が追加される行は、必ずしも矛盾ではなく「記述の粒度差」である。全18行は [`subunit_review.csv`](../runs/20260918_1117_chembl_acquisition/tables/subunit_review.csv) に記録した。

## 検証

- `MPLBACKEND=Agg .venv/bin/python -m pytest -q`: **20 passed**（既存4件＋追加16件）。macOS GUI backendによる最初のテストプロセス終了を避けるためheadless backendを明示した。
- ID必須、重複ID拒否、立体識別子/親ID不一致、subunit抽出、欠損/censoring、ページ不足・異常応答・外部host拒否、raw改変検知、失敗run明示、offline replayと元run不変性を検証。
- 本runのraw94応答とmanifestに記録した全artifactのSHA-256を検証。
- 本runからのoffline replayで、processedの全10ファイル（5 CSV＋5 JSONL）がbyte単位で一致。検証用replayはリポジトリ外の作業用ディレクトリで実施した。
- 既存の全追跡CSVについて、baseline inventoryのSHA-256と一致を確認。過去run・IFP・薬理・phenotype出力は未変更。今回既存ファイルで変わるのはREADMEとLATEST_RUNポインタのみ。

## 次の判断

各追加剤について、species、target/subunit、assay条件、endpoint、unit、relation、文献を揃えた比較可能性の確認が必要。今回は取得とQCを完了し、追加docking、phenotype自動label化、ML学習は実施していない。

[実装・再実行手順](DATA_ACQUISITION.md)と[既存解析レビュー](REPOSITORY_REVIEW.md)を参照。
