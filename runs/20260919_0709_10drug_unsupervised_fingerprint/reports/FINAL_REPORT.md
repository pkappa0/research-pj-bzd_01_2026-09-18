# 10-drug unsupervised fingerprint — preflight blocked

This is an audit/preregistration run, NOT a completed ten-drug analysis.

## Audit outcome / 不足項目
追加6剤Dockingは未実行。既存4剤と同一環境をrepositoryから再現できないため、依頼Step1のゲートで停止した。Vina/RDKit/PLIP/Open Babelの旧実行version、vendored PLIPの実効config、元Vina log/pose出力および全poseのPLIP reportが必要。requirements.txtは版固定なしでこれらを含まず、work/とdocking/はgitignore対象。現在Vina v1.2.7が存在することは旧versionの証明ではない。ローカルvenvにはRDKit/PLIP/Open Babelもない。別versionをインストールして同一pipelineと呼ぶことはしていない。

## 確認できた条件
6HUP ABCDE / 9CTJ CDE、保存済みPDBQT、box、seed20260917、exhaustiveness8、num_modes9、energy_range4、全返却pose採択を確認。boxは20 Å角で共通、中心はoutputs/docking_boxes.csvをそのまま使用する設計。ligandはhash seed→AddHs/EmbedMolecule→MMFF（例外時UFF）→RemoveHs→Gasteiger簡易writer→TORSDOF0。受容体電荷は0、pH指定なし。EmbedMoleculeの未指定defaultがversion依存であるため、コメントのETKDGだけから実効設定を推測しない。新薬のために物理的に改善した前処理へ置き換えると旧4剤との条件差になる。

既存4剤100 candidate feature行は元contact表から頻度を再計数して完全一致、72 pose行（各drug/block 9）を確認。ただしこれは保存CSV間の整合性であり、存在しないPLIP元reportによる全pose成功監査の代わりではない。未観測を無条件0とした10剤表は作成していない。

## 化合物構造
10剤のChEMBL ID、SMILES、InChIKey、取得元hashを既存取得runから抽出。追加6剤も別sourceへ切替不要。lorazepam / temazepam / zopicloneの記録SMILESには立体指定がなく、単一生成conformerをracemateの代表と断定できない。未指定立体・microstateの扱いは旧環境回収後も明示する必要がある。

## Preregistration
config/unsupervised_analysis_config.jsonは新規Docking/PCA/clustering前に保存。Bは各残基のtypeをcollapseしたpose union frequency、Cはtype別頻度＋P/T、EはBのα1/α2連結。A binary、D geometry付きも事前定義。primaryはEuclidean＋average linkage、k=3、center-only PCA。secondaryのcosineとk=2/4も事前指定して全結果を報告し、外部annotationに合う条件を選ばない。PCA/clusteringは各表の全10剤共通観測列のみ、pairwise similarityは共通feature数も保存。0補完は禁止。heatmapはtype別C行とprimary B/E dendrogramを表示する計画。geometryには別単位と欠損maskを保持する。

ここでは薬理activityやmouse phenotypeをロードしていない。構造結果hashをfreezeしてからannotation処理を行うゲートを設定。既存会話で知られている薬理知識から解析条件を最適化しない。単一seedのmetric感度を独立Docking再現性の証拠とはしない。

## A. Structural pattern
**indeterminate / 未評価**。追加6剤のIFPがなく、10剤cluster/PCAは実行できない。4剤の一致監査を10剤の結果へ一般化しない。

## B. External correspondence
**indeterminate / 未評価**。構造結果未固定のため、外部annotationの統合は未実行。

## C. ML expansion rationale
**indeterminate**。現段階の不足は再現環境/provenanceの問題であり、科学的仮説を支持も反証もしない。10剤matrix、similarity、cluster assignments、PCA図は未生成。欠損6剤を0で埋めた図やダミー図は作らない。

## 再開に必要な証拠
旧実行環境のpackage lock/conda export/container digest、Vina --version記録、work/plip_vendorと実効config、元strict3/strict4 docking job directories（Vina log、poses.pdbqt、per-pose PLIP XML）。これらが回収不能なら、完全同条件の追加解析という条件を満たせない。条件変更を伴う別設計は、ユーザーによる研究方針の変更を受けてから別runで扱う。
