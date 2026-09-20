# C1 direction consistency — median-ordered descriptive visualization

既存forest plotと同じ数値を使用した別図。利用者の今回の明示的な指定により、`median_rho_across_c1_terms`の降順で48 featureを並べる。**並べ替えはdescriptive visualizationのみで、feature selectionではない。** 同値は元のfrozen feature順を維持し、NAがある場合は末尾へ置く。元のfrozen-order figureは上書きしない。

## Inputs and unchanged values

Source run: `runs/20260920_c1_direction_consistency/` (commit `3429d776fe77c41a43f2f1d0dcaee7efccdfb862`).

- `data/c1_direction_consistency_summary.csv`: original medians, min/max, estimable-term and sign counts.
- `data/feature_display_key.csv`: exact feature labels.

上記CSVのbyte-identical snapshotを本runの`data/`へ保存。`input_manifest.json`に元pathとSHA256を記録。元のANALYSIS_CONTRACT.mdもsnapshotする。既存相関、C1 summary、PLIF、Clinical fingerprintを再計算・変更せず、表示順だけを変更する。受容体を跨いだ平均やcomposite scoreは作らない。

## Figure

点は元のmedian rho、横区間は元のmin–max rho。横軸は−1〜+1に固定し、0の基準線を置く。行ラベルにF IDとα1/α2を明記し、青＝α1、橙＝α2として所属を補助表示する。全48 featureを保持。

右注記は今回指定された3分類を使用する：4語すべて正なら`4/4 positive`、4語すべて負なら`4/4 negative`、それ以外は`mixed`。前図の3/4表示は今回mixedへまとめる。F29/F48のように4語未満しかestimableでない場合は`mixed (2/4 estimable)`などとcoverageを明記し、欠測をゼロや4/4一致とみなさない。全語未定義なら`NA (0/4 estimable)`。入力のsign countsは変更しない。

区間はclinical terms間のrhoのvariationであり、confidence intervalではない。統計単位はdrug。Association/correspondence only; no causal or mechanistic inference. 順位は統計的有意性やfeature重要度を表さない。cutoff・選別・hit callingを実施しない。

## Reproduce and output paths

Python + numpy + pandas + matplotlib:

```bash
python runs/20260920_c1_direction_consistency_ranked/scripts/plot_c1_direction_consistency_ranked.py
```

- `figures/c1_plif_direction_consistency_ranked.png`
- `figures/c1_plif_direction_consistency_ranked.svg`
- `data/figure_display_annotations.csv`: sorted display metadata with unchanged numeric values.
- `QC_REPORT.json`: source hashes, row order and checks.

PNG 3240×5040 (240 dpi); SVG preserves vector points/intervals and text. Both are separate files in a new run; no historical outputs are overwritten.
