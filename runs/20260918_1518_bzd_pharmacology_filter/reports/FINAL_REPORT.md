# BZD-site pharmacology filter — 20260918_1518_bzd_pharmacology_filter

入力: `runs/20260918_1117_chembl_acquisition/processed/activities.jsonl`（ChEMBL 37、434 activity、10剤）。
取得済みアーカイブだけを使用。追加API取得・ML・Docking・phenotype予測は実施していない。

## 結果

| category | activity数 |
|---|---:|
| PRIMARY | 70 |
| SECONDARY | 193 |
| REVIEW | 157 |
| EXCLUDE | 14 |
| 合計 | 434 |

434行すべてを1カテゴリへ分類し、元activity IDと1対1対応。PRIMARYは直接比較の**候補**であり、独立した検証済み教師ラベルではない。
PRIMARYには数値を保持したcensored値、potential_duplicate、assay organism欠損も含む。厳しいpoint比較条件を満たす件数を別に示す。

- PRIMARY subtype: `{"alpha1":41,"alpha2":29}`。alpha1_alpha2 / other_explicit_subtypeはPRIMARY 0件。
- 全件endpoint: `{"binding":286,"functional_efficacy":76,"functional_potency":14,"other":58}`
- PRIMARY endpoint: `{"binding":58,"functional_efficacy":11,"functional_potency":1}`
- 全件species: `{"bovine":16,"human":101,"mouse":5,"rat":159,"unknown":153}`
- PRIMARY species: `{"human":29,"rat":21,"unknown":20}`
- PRIMARY point比較候補: 27件、comparison group: 52件。

## Compound coverage

usableはPRIMARYかつ有限値・relation `=`・報告単位あり・assay species既知・potential_duplicateなし。
各α1/α2 × binding/potency/efficacyの件数、group数、human有無、review有無は `compound_pharmacology_coverage.csv` に保存。

| compound | ChEMBL ID | IFP済 | PRIMARY | α1 usable | α2 usable | binding | functional | human | REVIEW行 | flags延べ数 | 記録条件pair数 |
|---|---|---|---:|---:|---:|---|---|---|---:|---:|---:|
| diazepam | CHEMBL12 | True | 34 | 12 | 5 | True | True | True | 90 | 151 | 3 |
| lorazepam | CHEMBL580 | False | 0 | 0 | 0 | False | False | False | 2 | 3 | 0 |
| clonazepam | CHEMBL452 | False | 0 | 0 | 0 | False | False | False | 5 | 7 | 0 |
| alprazolam | CHEMBL661 | True | 2 | 1 | 1 | True | False | True | 11 | 28 | 1 |
| midazolam | CHEMBL655 | False | 0 | 0 | 0 | False | False | False | 2 | 3 | 0 |
| temazepam | CHEMBL967 | False | 0 | 0 | 0 | False | False | False | 3 | 4 | 0 |
| triazolam | CHEMBL646 | True | 2 | 0 | 0 | False | False | False | 2 | 3 | 0 |
| zolpidem | CHEMBL911 | True | 32 | 3 | 5 | True | True | True | 33 | 69 | 0 |
| zopiclone | CHEMBL135400 | False | 0 | 0 | 0 | False | False | False | 3 | 5 | 0 |
| zaleplon | CHEMBL1521 | False | 0 | 0 | 0 | False | False | False | 6 | 8 | 0 |

α1/α2双方にPRIMARYあり: **diazepam, alprazolam, triazolam, zolpidem**。
α1/α2双方にpoint比較候補あり: **diazepam, alprazolam, zolpidem**。
同一compound・文献・species・β/γ・endpoint・単位・記録条件のα1/α2 pair候補あり: **diazepam, alprazolam**。

pair候補はassay IDをまたぐため、α1/α2表記だけをマスクしたdescriptionと、保存されたcell/tissue/strain/parameter等の完全一致を要求する。
未記載条件の一致を証明するものではない。descriptionの表記差で見逃す可能性があり、pairなしは生物学的比較不能の証明ではない。
comparison group自体はassay ID・文献・species・full composition・type・unit・relationごと。compound名はkeyに含めず、平均・代表値・比は計算しない。

## 既存IFP 4剤と追加6剤

diazepam・alprazolam・triazolam・zolpidemは既存IFP群として保持。triazolamのPRIMARYはassay organism欠損のためpoint比較候補に数えない。
zolpidemのα1/α2 usable行が存在しても、厳密な記録条件pairが成立するとは限らない。
alprazolamのpair候補はβ2背景を含み、既存IFPのβ3背景と同一構成とは扱わない。pharmacology内のpairとIFP構造との一致は別途確認が必要。
lorazepam・clonazepam・midazolam・temazepam・zopiclone・zaleplonは未IFP群。骨格・compound名による除外は行っていない。
現在のルールで追加Dockingを検討するための記録条件pairを持つ未解析compound: **なし**。
該当しない薬剤に順位や科学的best drugを付与しない。`additional_docking_coverage_rank` の欠損は評価不十分を意味する。
zaleplonは原著確認の候補。current modulationのEC50なのにassay_type BとされるためREVIEW。またα1β2γ2とα2β3γ2はβ背景が異なり、注釈を確認してもそのままmatched pairにはできない。
zopicloneのnative receptor/current efficacyも構成とassay分類の確認が必要。残る4剤はsubtype-completeな比較データの追加探索が先となる。

## Major review issues

| review flag | activity数（重複計上あり） |
|---|---:|
| assay_type_description_conflict | 82 |
| bzd_site_relevance_unconfirmed | 55 |
| target_assay_organism_conflict | 29 |
| concentration_endpoint_context_ambiguous | 18 |
| assay_origin_or_measurement_uncertain | 17 |
| target_description_organism_conflict | 13 |
| repeated_subunit_sequence_or_concatemer_review | 11 |
| variant_construct_requires_review | 10 |
| missing_endpoint_unit | 9 |
| missing_standard_value | 6 |
| chembl_data_validity_comment | 5 |
| missing_or_unrecognised_relation | 5 |
| endpoint_interpretation_ambiguous | 4 |
| target_assay_alpha_conflict | 3 |
| target_component_assay_alpha_conflict | 3 |
| binding_functional_description_ambiguous | 2 |
| target_assay_beta_conflict | 2 |
| target_component_assay_beta_conflict | 2 |
| affinity_endpoint_without_binding_description | 2 |
| log_endpoint_scale_requires_review | 1 |
| concentration_endpoint_vs_fixed_percent_description | 1 |
| compound_relative_test_concentration | 1 |

α1対α4、α1対α5、β1対β2の矛盾は一般的なfamily別集合比較で検出。activity IDに依存した分類分岐はない。
target/assay species不一致、host organism由来の可能性、unknown originの記述も自動修正しない。
assay organism欠損は全件中124件。target/proseから補完せずunknownとする。矛盾時も解析用speciesはunknown、各元speciesは保持。
B/F assay分類は粗い注釈なので、不一致は実験が誤りという断定ではなく原著確認要求。REVIEW優先のため、off-siteと判断できる行でも別の矛盾があればREVIEWに残る。
component metadataのGROUP MEMBERは受容体の共集合構成と解釈しない。single-subunit targetやgeneric/nativeはPRIMARYに昇格させない。
functional_efficacyは広いresponse分類。固定濃度応答と真のEmaxはendpoint_detailで区別し、同じgroupに入れない。pChEMBLへの統一やKi/IC50/EC50/Emaxの数値統合は行わない。

## QC・provenance・再現

QC: **PASS**。18チェックを通過。入力raw gzipと元manifestのartifact hashも検証。
tests: `{"errors":0,"failures":0,"skipped":0,"source_sha256":"26c299aafe6d7167a55abb2607716d6dec3d4446d146635a62d4408c6e0c1957","status":"PASSED","tests":70}`。詳細は `logs/pytest.xml`。
`reports/QC_REPORT.json` はID一意性、元全fieldの完全保持、null・relation・censoring・species・composition・group整合性を記録。
`raw/lineage/` は入力processed JSONL、registry、元manifest、raw request ledgerのbyte-identical snapshot。
元APIレスポンスはsource runの `raw/responses/` に保持し、コピーし直さない。元source_file / JSON pointerは **source run基準**。
`processed/bzd_activity_classification.jsonl` はnullと元数値文字列を保持する正本。CSVの空セルは欠損（0ではない）、list/dictはJSON表現。
各activityへsource run、元processed行番号・ファイルhash、rule versionを追記。元target/component/assay/document情報とURL・取得時刻・raw hashは失わない。
実行時コード、config、tests、依存関数を `config/code_snapshot/` に保存。manifestのgit_commitは実行前HEADであり、実行コードの厳密な識別はsnapshot SHA-256による。
既存tracked filesは実行前後のhash一致を確認。旧run、IFP、phenotype、`runs/LATEST_RUN.txt` を更新しない。

再実行: `MPLBACKEND=Agg .venv/bin/python -m pytest -q --junitxml=/tmp/bzd-tests.xml` の後、
` .venv/bin/python -m src.bzd_pharmacology_filter --config config/bzd_pharmacology_filter.json --test-results /tmp/bzd-tests.xml`。
再実行は別runを生成する。ルールは保守的な機械抽出であり、negationや複雑なconstructを完全には解釈できない。PRIMARYも原著とIFP構造のspecies/subunit/条件を照合してから使用する。

## 判定設計の参照

- [ChEMBL data FAQ](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions): assay分類、target注釈、potential_duplicateを別々の情報として扱う。
- [BZD-site構造研究](https://www.nature.com/articles/s41586-018-0255-3)、[GABAA構造薬理](https://www.nature.com/articles/s41586-018-0832-5): α/γ界面と明示構成を重視。
- [TBPS binding study](https://pubmed.ncbi.nlm.nih.gov/3035434/): channel-site readoutをBZD-site affinityと混同しない。
