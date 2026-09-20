# C1 correspondence arranged by residue properties

既存forest plotと同じmedian/min/max rhoを、残基性質の表示グループごとに配置した別図。全48 featureを保持し、各グループ内は直前のranked図と同じmedian rho降順、同値は元のfrozen順。**表示上の並べ替えのみでfeature selectionではない。** 相関・C1 summary・frozen PLIF・Clinical fingerprintの変更や再計算を行わない。受容体や残基群を跨ぐ平均・合成scoreも作らない。

## Display groups

| 順序 | 表示分類 | 本図で観測される残基 | feature数 |
|---|---|---|---:|
| 1 | 芳香族 / Aromatic | PHE, TYR, TRP | 11 |
| 2 | 非芳香族疎水性 / Hydrophobic, non-aromatic | ALA, VAL, ILE, PRO | 7 |
| 3 | 親水性・極性非荷電 / Hydrophilic, polar uncharged | SER, THR, ASN, GLN | 15 |
| 4 | 酸性・解離可能 / Acidic, ionizable | ASP, GLU | 9 |
| 5 | 塩基性・解離可能 / Basic, ionizable | LYS, ARG, HIS | 6 |

feature数は薬剤数ではない。同じ残基種類・部位に複数receptorのfeatureが存在しても、それぞれ元featureとして保持する。α1は青、α2は橙、各行にもα1/α2を明記。

性質は本来重なるため、これは互いに排他的な**表示用分類**であり、一義的な物理化学的分割ではない。酸性・塩基性の側鎖も広い意味では親水的だが、電荷に関わる性質を区別するため独立群にした。TYRはフェノール性OHを持つ極性のある芳香族残基で、ここでは芳香族へ置く。HISのイミダゾール環も芳香族だが、表示上は塩基性・解離可能へ置き、`aromatic_sidechain_flag=true`を別欄に残す。TRPのindole nitrogenも備考として記録する。実際の受容体環境でのprotonation/電荷、溶媒露出、hydropathy値を計算・推定したわけではない。

`residue_class_rules.json`は表示用対応規則、`data/residue_class_mapping.csv`は全48 featureの具体的分類、芳香族flag、補足。規則中に本図で未観測の標準残基があっても、人工的なfeature行は追加しない。

## Sources and outputs

数値入力は既存fixed-order forest run `runs/20260920_c1_direction_consistency/` の:

- `data/c1_direction_consistency_summary.csv`
- `data/feature_display_key.csv`

byte-identicalなコピーを本runの`data/`に保持し、元pathとSHA256を`input_manifest.json`に記録。元のANALYSIS_CONTRACT.mdも保存。本runの追加は表示分類と行順のみで、解析の比較軸・feature定義は変わらない。

Outputs:

- `figures/c1_plif_direction_consistency_residue_class.png`
- `figures/c1_plif_direction_consistency_residue_class.svg`
- `data/figure_display_annotations.csv`: 元数値と並べ替え後のラベル/分類/方向注記。
- `data/residue_class_mapping.csv`, `residue_class_rules.json`, `QC_REPORT.json`.

PNG 3240×5520 (240 dpi); SVGはvector/textを保持。従来のfrozen-order図とmedian降順図は上書きしない。

## Reading the figure

点は既存median rho、区間は既存min〜max rho、横軸は−1〜+1。右注記は4/4 positive、4/4 negative、それ以外mixed。4語未満はestimable数を明記。欠測の補完はしない。区間はC1 clinical terms間の変動であり信頼区間ではない。Observational unit = drug. Association/correspondence only; no causal or mechanistic inference. 分類ごとの並びから因果残基、薬効機序、統計的有意性を主張しない。

## Reproduce

Python + numpy + pandas + matplotlib:

```bash
python runs/20260920_c1_direction_consistency_residue_class/scripts/plot_c1_direction_consistency_residue_class.py
```

分類・48-row完全性・群内降順・元median/min/max一致・入力hashを検証。順位付けのための再相関計算や閾値設定はない。
