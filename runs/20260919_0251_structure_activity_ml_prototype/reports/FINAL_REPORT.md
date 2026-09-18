# Structure–activity ML dataset prototype

## Research question
Residue interaction fingerprintに、公開DBのpharmacological activity profileを説明・予測する情報が含まれるか。主解析は **same receptor → different drugs**。本runはX–y形式の構築とQCのみで、学習・CV・feature importanceを実行しない。

## Current evidence
薬剤2行、構造Xはα1 27列＋α2 27列（36 frequency/type、18 geometry）。同一構造block内で薬剤差があり、同一記録assay context内でもKiの数値差がある。対応表を構築できたことは予測性能の証拠ではない。

| Receptor-fixed pharmacology | Diazepam Ki (nM) | Alprazolam Ki (nM) | Δ(alp−dia), nM | Dia/alp | Assay |
|---|---:|---:|---:|---:|---|
| alpha1_beta2_gamma2 | 14.0 | 0.8 | -13.2 | 17.5 | CHEMBL5146040 |
| alpha2_beta2_gamma2 | 20.0 | 0.6 | -19.4 | 33.3333 | CHEMBL5146041 |

Kiは小さいほど報告された結合affinityが高い。すべてhuman、relation =、nM、同一document CHEMBL5143601。ただし不確実性推定はなく、有意差・選択性ラベルへ変換しない。activity IDs、assay descriptions、raw source pointersはtables/pharmacology_targets.csvとprovenance/selected_activity_full.jsonlに保持。

### Mechanically selected presentation features
| Feature | Diazepam | Alprazolam | Δ(alp−dia) |
|---|---:|---:|---:|
| α1 γ2 TYR58 hyd | 0.777778 | 1.000000 | +0.222222 |
| α1 γ2 TYR58 π | 0.444444 | 0.444444 | +0.000000 |
| α1 γ2 TYR58 π/P | 0.444444 | 0.444444 | +0.000000 |
| α1 γ2 PHE77 hyd | 0.777778 | 0.777778 | +0.000000 |
| α1 γ2 PHE77 π | 0.222222 | 0.444444 | +0.222222 |
| α2 γ2 TYR58 hyd | 0.333333 | 0.444444 | +0.111111 |
| α2 γ2 PHE77 hyd | 0.666667 | 0.666667 | +0.000000 |
| α2 α LYS156 hyd | 0.222222 | 0.555556 | +0.333333 |
| α2 γ2 ASN60 hyd | 0.333333 | 0.444444 | +0.111111 |

α1ではTYR58 hyd、PHE77 πに増加があり、TYR58 π/PやPHE77 hydなどは薬剤間で同頻度でも表示対象となる。α2ではLYS156 hydの増加（2/9→5/9）、TYR58 hydとASN60 hydの増加を保持。これは受容体間比較ではない。残存頻度が高いことと薬剤識別力があることは別である。

## Raw / reduced policy
Rawには定義した54列をすべて保持。完全な既存IFP行（geometryのmin/max/median、観測数、pose/source IDsも含む）はprovenance/source_ifp.jsonlに保存。図とrawのgeometry列はcontact-row平均centroid distance / angle / offsetで、観測contact条件付き・pose独立ではない。P/T混合時の単一平均には限界があり、型別頻度と元contact provenanceを併用する。
縮約は両薬剤frequency <0.40かつ絶対差≤2/9を除外する表示ルール。整数contacting-pose数/pose分母のFractionで境界を厳密評価する。9 frequency列と、残存π頻度に紐づく6 geometry列を保持する。Geometryにfrequency閾値を適用しない。その他39列はrawに残る。欠損頻度は除外せず保持。これは統計的有意差でもML用feature selectionでもない。

## X–y context and missingness
構造α1はα1β3γ2L full pentamer、α2は9CTJ由来provisional β3/α2/γ2 local C/D/E。薬理はα1/α2β2γ2であり、観測された4つのX–y対応はすべてcontext_mismatch。薬剤間assayのexact_context_matchとは別のQC軸である。未観測yはnot_assessed_missing_activityとし、存在しない測定にexact/partial/mismatchを捏造しない。
wide datasetはraw54 X列＋観測Ki 2列＋未取得endpoint 12列（各blockのKd、other binding、EC50、IC50、Emax、response）を保持。未取得はNA、unit/relation/assayもNAであり、ここでのreceptor/speciesは将来取得の要求context。取得済みdrug片側のみの活動をこの対応表のyへ転用しない。該当15 eligible unpaired recordsも別provenanceに保持する。IC50はassay機序でbindingの場合もあるため、将来もendpoint名だけでfunctional potencyと分類しない。
0は既存監査済みcontactなし。NAは未観測/未定義であり0補完しない。ml_dataset_prototype.csvのqc列・drug ID・ChEMBL IDは予測featureではない。feature_schema/target_schemaとmanifestがX/y列の契約である。CSVではNAを明示、missingness_mask.csvも保持。

## Future model comparison manifest
Model0: 各blockの9 archived Vina poseの最小score。Model1: 6対象siteで指定interaction typesのいずれかがあるbinary。Model2: 同じsiteのcontact pose集合の和集合/9（type頻度の和ではない）。Model3: type別frequencyとP/T。Model4: Model3＋geometry。Model5: α1/α2のModel4を順序固定で連結。Model0–4はblock別に定義し、Model5の増分は事前指定single-block baselineと比較する。baseline_features.csvとfeature_set_manifest.csvに実列名を保存。Model1/2は対象6site・選択interaction types内のunionであり全残基網羅ではない。
すべて同じy endpoint、assay/context、評価drug集合・splitで比較する。欠損geometry等の処理は将来training fold内のみで決める。Reduced表示ルールを学習feature選択に流用しない。endpointを平均/統合しない。

## Not claimed
- 特定残基がKi差を原因として生じさせる。
- 現時点の2剤で予測性能がある。
- すべての薬理差を構造Fingerprintだけで説明できる。
- 小さいinteraction差すべてに生物学的意味がある。
同方向に動くfeatureがあっても、独立drug n=2、単一pose生成seed、同一化学系列、構造/薬理context不一致により対応情報の有無は判断不能。9 poseを独立drugサンプルと扱わない。

## Next step
最優先は同一組成・species・assay機序/endpointで複数追加薬剤のyを確保し、対応する構造Xを作ること。β2対応構造またはβ3比較薬理でcontext差を分離する。十分なdrug数と独立評価集合を確保した後、Model0〜5を同じendpointで比較する。反証可能な問いは「同一contextの未知drugで、type/geometry追加やmulti-receptor化がdocking score/単純contact頻度よりheld-out予測誤差を改善するか」。改善がなければ追加表現の情報価値仮説を支持しない。現時点では拡張の実行可能性は支持、予測上の価値は判断不能。

## Provenance / reproduction / QC
ANALYSIS_CONTRACT.mdのsnapshot/hashを保存し、既存runは変更しない。source_inventory.json、selected_activity_full.jsonl、source_ifp.jsonlから元runとChEMBL取得レスポンスへ遡れる。
`MPLBACKEND=Agg .venv/bin/python -m pytest tests --junitxml=<new-path>/tests.xml`
`MPLBACKEND=Agg .venv/bin/python -m src.structure_activity_ml_prototype --tests-xml <new-path>/tests.xml`
テストとQC結果はreports/tests.xml、reports/QC_REPORT.json。図はPNG/PDF、numeric scaleをfrequency・geometry各単位・Kiに分離する。既存runと値を照合し、block非平均・missing保持・endpoint分離・context flags・source IDs・学習未実行を確認。
