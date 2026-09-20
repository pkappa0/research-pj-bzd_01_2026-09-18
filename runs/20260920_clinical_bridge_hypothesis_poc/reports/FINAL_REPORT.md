# 三層 working hypothesis 予備検証 — 2026-09-20

## 結論

固定済み48 structural featureと、独立に事前定義した6軸14語の臨床fingerprintを、同じ10剤identityで接続した。教師ありML・新規docking・feature再選択は実施していない。今回の臨床指標は **all-role openFDA報告文書に基づく未調整log(ROR)** であり、primary suspect解析ではない。

- **A. Clinical fingerprint feasibility: supported（取得方法の範囲内）**。10剤中10剤で安定性基準を満たす値を作成。6軸を保持し、14語を単一scoreへ統合しない。用語ごとの低countは下表参照。厳密なprimary-suspect cohort、完全な個票共変量・重複症例監査まで実現したという意味ではない。
- **B. X ↔ Y correspondence: indeterminate**。45 drug pairsの記述的Spearman rho=-0.2000。9999回の薬剤ラベル置換による両側tail fraction=0.4548。45組を独立n=45とする検定はしていない。全estimable値を含める感度解析rho=-0.1593、1剤除外rho範囲=-0.3578〜-0.0770。
- **C. Intermediate phenotype data-linkage feasibility: supported**。既存in vivo evidenceを同一drug identityで接続可能。比較可能な数値anchorはNishino2008のdiazepam/triazolamの2剤のみ。中間機能そのものの支持ではない。
- **D. Three-layer working hypothesis: indeterminate**。表現とデータ連結の実現可能性については根拠があるが、全体のX↔Yは正の対応を示さず、M/Y anchorも一貫しなかった。**in vivoが構造と臨床を接続する中間表現として機能するという科学的主張はindeterminate**。相関、mechanism、mediation、因果経路の証明ではない。

MLへ直ちに進む根拠は不足。まず薬剤数・比較可能なM・臨床集計の妥当性を増やし、固定仮説の反証可能な追試を優先する。

## 再現性と解析契約

run: `runs/20260920_clinical_bridge_hypothesis_poc/`

source repository commit: `cff06d2`（完全hashはstructural manifest）。成果commitは配布物のCOMMIT.txt / Git履歴を参照。commit自身のhashをそのcommit内へ埋め込む循環は避ける。

元のANALYSIS_CONTRACT.mdをsnapshotし、利用者が指定したX↔Yへの科学的比較軸拡張をconfig/CONTRACT_ADDENDUM.mdに明記。α1/α2 blockを保持し、各blockは反復標本ではない。構造入力は20260919_0724_10drug_standardized_redockingの48列をそのままコピーしSHA256・列順・元commitを保存。α2 provisional construct、β3構造と他薬理系の不一致に関する元の制約は残る。

STRUCTURAL_RESULTS_FROZEN.json → clinical config / raw acquisition → Y-only QC → CLINICAL_RESULTS_FROZEN.json → 対応解析の順。analysis scriptはX/Y hash一致を必須にする。Y結果によるterm選び直しはない。

## Clinical sourceと母集団

[openFDA drug event API](https://open.fda.gov/apis/drug/event/)のFAERS報告文書。`receivedate:[20040101 TO 20251231]`を全薬剤・背景で共用。背景は20,018,532文書。取得日時、完全query、HTTP status、レスポンスhashをrawに保存。APIの更新日時は各response metaに保存（本取得の背景responseは2026-07-30更新）。2026年の部分年は含めない。

薬剤同定は`openfda.generic_name.exact`または`medicinalproduct.exact`のgeneric表記をOR結合。ブランド別countを足さない。ChEMBLのcompound identity / InChIKeyは固定runにリンクし、API由来のbrand / other synonymを別表へ保存。非標準ブランド・塩・複合剤など、exact identityで拾えない報告を無理に採用しない。zopicloneとeszopicloneは統合しない。この方法は特異性寄りで、全処方・全報告の完全回収を保証しない。generic-only総数は同定coverage QCとして保存し、別のsuspect解析とは呼ばない。

openFDAの[公式field definition](https://open.fda.gov/fields/drugevent.yaml)ではdrugcharacterizationはsuspect / concomitant / interactingで、PS/SSは分離されない。また報告内に複数薬剤があるため、薬剤名条件とrole条件をreport-level ANDにしただけではその薬剤自身のroleを保証できない。従ってprimary-suspectとall-suspect感度解析は未実施。all-roleを代替primaryとしたことを明示する。

同定されたreport documentを数え、同一report中の複数表記をORで重複加算しない。公式field definitionではAPIは最新report versionを返す。一方、別IDとして存在する重複症例を個票で独立除去したわけではない。

性別・年齢・国・受領日・seriousnessの薬剤別周辺分布をclinical_demographic_marginals.csvへ保存。上位100値の集計なので分布が切れる場合がある。ageはunitとjoint集計できておらず年齢yearsへ変換しない。薬剤arrayのrole/brand周辺値には併用薬の値も含み、その薬剤への帰属を主張しない。**drug×eventごとの個票共変量、完全なquarter cohort、調整RORは未作成**。clinical_raw_counts.csvのALL表記は全層集計であり、取得した個人属性を意味しない。生responseの例示recordは完全個票cohortではない。

## 用語事前定義とY-only QC

| Axis | 1語以上安定値がある薬剤 | 採用term |
|---|---:|---|
| C1 | 10/10 |Ataxia / Balance disorder / Coordination abnormal / Gait disturbance |
| C2 | 10/10 |Fall |
| C3 | 10/10 |Somnolence / Sedation |
| C4 | 10/10 |Dizziness / Vertigo |
| C5 | 10/10 |Muscular weakness / Hypotonia |
| C6 | 10/10 |Amnesia / Memory impairment / Cognitive disorder |

openFDA reactionmeddraptの正確なsource文字列を個別featureとする。対応の確定できないabnormal gait、tendency to fall、injury secondary to fall、drowsiness、muscle relaxationは未mapping行として残し、勝手な同義語unionはしない。VertigoはDizzinessと分離。これはデータソース文字列へのmappingであり、licensed MedDRA ontologyのPT/LLT照合を完了したとの主張ではない。採用14語は背景APIで存在を確認できた。

| Drug | All-role documents | 安定feature数 | a<5のセル |
|---|---:|---:|---:|
| diazepam | 118,049 | 14/14 | 0 |
| alprazolam | 200,932 | 14/14 | 0 |
| triazolam | 3,893 | 12/14 | 2 |
| zolpidem | 49,105 | 14/14 | 0 |
| lorazepam | 174,733 | 14/14 | 0 |
| clonazepam | 155,429 | 14/14 | 0 |
| midazolam | 24,085 | 14/14 | 0 |
| temazepam | 36,961 | 14/14 | 0 |
| zopiclone | 42,339 | 14/14 | 0 |
| zaleplon | 1,544 | 10/14 | 4 |

取得失敗は例外として停止し、zero countへ変換しない。APIの明示NOT_FOUNDのみ0件。背景存在または薬剤総数がない場合はRORをNA。a<5はunstableとして主YではNA、raw RORには残す。低countを除外することで薬剤ごとの表現の精度差が生じる点にも注意する。

## RORと対応解析

a=drug+event、b=drug+not event、c=background event−a、d=background total−a−b−c。drug群とnot drug群を同じ背景から構成する。複数target drugを含む報告はそれぞれの薬剤の解析に含まれ、薬剤間の統計的独立性も保証されない。

ROR=ad/bc、logROR=log(a)+log(d)−log(b)−log(c)。ゼロcellがあれば全4cellへ0.5加算。log scale SE=sqrt(1/a+1/b+1/c+1/d)、nominal 95% CI=exp(logROR±1.96SE)。CIは報告バイアス、症例重複、交絡を含む総不確実性ではない。

Xは48 frequency列のcosine、Yは固定14語のlogRORのcosine（pairwise complete、7語以上必要）。欠測をゼロ補完しない。45組それぞれの共有feature数を保存。cosineはprofile方向の近さであり、同一の副作用頻度やabsolute riskの近さを意味しない。全般的報告増加とaxis特異的な構造対応も分離できない。

feature×termは全48×14=672組を保存し、n_effective、missingness、rhoを記録。n<5はexploratory-insufficient、constant featureはrho=NA。p-value rankingを作らない。

**最大の記述的候補**: `alpha2|BZD_GAMMA2_192|ASP|any_contact_frequency` × `Muscular weakness`、rho=0.8964、n_effective=10。これは672比較中の事後的な最大絶対値の記述であり、独立な発見・残基機序・再現性の保証ではない。同率候補は全表に残る。共通交絡や構造feature間相関によっても生じ得る。

## Diazepam / triazolam三層例

Xは全featureのtriazolam−diazepam差を表示。Mは同一研究、同一経口dose 2/5 mg/kg、同一time 15/30/60/90分の8context。Triazolamのrotarod failure fractionは7条件で高く、1条件で同値。n=2薬剤であり、8条件を8薬剤対とは扱わない。Mの未報告値はNA。

Yの差は以下。正はtriazolamのlogRORが高いことを表す。

| Clinical term | logROR triazolam−diazepam |
|---|---:|
| Ataxia | -0.057 |
| Balance disorder | -0.447 |
| Coordination abnormal | NA (a<5) |
| Gait disturbance | -0.206 |
| Fall | +0.005 |
| Somnolence | -0.385 |
| Sedation | -1.059 |

この方向比較は未調整の報告傾向と動物assayの記述を並べただけである。今回、Ataxia、Balance disorder、Gait disturbance、Somnolence、SedationはいずれもtriazolamのlogRORが低く、rotarod障害が高い方向とは逆だった。Fall差は約0.005でほぼ同値、Coordination abnormalは低countでNA。従って一貫したmotor/clinical方向一致は観察されない。Xには事前に検証された「motor impairment方向」のscoreがないため、構造差の正負だけから三層の機序的方向一致は判定できない。結論は **data linkageは可能、三層の機能的整合性はindeterminate**。n=2のcorrelation / mediation / causal pathは計算しない。

## MLへ進む前に不足するもの

1. より多い独立薬剤と、固定したX/Y定義を使う外部検証。化学構造の類似性など簡単な代替説明との比較。
2. FAERS bulk個票によるPS/SS・対象薬剤roleの正確な紐づけ、case重複監査、salt/brand coverage監査、同一期間/国/適応/併用薬の感度解析。
3. MedDRA versionとPT/LLT mappingの監査、および不均一な少数報告の不確実性評価。labelは別evidenceとして扱いFAERSと数値混合しない（今回label frequencyは使用していない）。
4. 同一assay/time/routeでdose・脳内曝露・species contextを揃えたMの多剤データ。P2/P3/P4を無理に補完しない。
5. 今回の最大相関をそのままfeature選択に使わない事前登録、holdout評価、交絡検討。n=10で教師ありMLへ移らない。

**最大のlimitation**は、構造・動物・自発報告が異なる測定過程であり、drug identity以外の臨床条件・曝露・適応・併用薬を揃えられていないこと。これにn=10、Mの比較可能な数値n=2、FAERSの報告偏りと共変量未調整が重なる。ROR≠risk ratio、報告frequency≠incidence、自発報告≠causal adverse effect、structural correlation≠mechanism、in vivo alignment≠mediation。

## 検証と成果物

code/validate.pyはfreeze hash、contingency再構成、ROR/CI、missingness規則、45 unique pairs、672 feature×term行、raw response hashを確認する。図5種をPNG/PDFで作成し、PDFを再renderして目視確認する。旧runと元入力は変更しない。結果の正負に合わせた条件変更は行わない。

- figures/three_layer_hypothesis.png / .pdf
- figures/clinical_phenotype_fingerprint_heatmap.png / .pdf
- figures/structural_vs_clinical_similarity.png / .pdf
- figures/structural_feature_clinical_phenotype_heatmap.png / .pdf
- figures/diazepam_triazolam_three_layer_example.png / .pdf

## 出典

- [openFDA FAERS overview](https://open.fda.gov/apis/drug/event/) — データの性質・非因果性・incidence推定不可。
- [openFDA field schema](https://open.fda.gov/fields/drugevent.yaml) — roleと最新report version。rawへsnapshot保存。
- [openFDA query parameters](https://open.fda.gov/apis/query-parameters/) — report search / count。
- [ChEMBL API](https://www.ebi.ac.uk/chembl/api/data/docs) — compound synonyms。各実queryはdrug_name_mapping.csv、JSONはraw。
- [Nishino2008](https://doi.org/10.1254/jphs.08107FP) — 既存runで監査済みTable1。今回新規数値抽出はしていない。
