# Frozen X/Y correspondence visualization

既存frozen X/Yをdrug identityで揃えた記述的可視化。**association/correspondence only; not causation or mediation**。

## Figure

`figures/frozen_xy_paired_heatmap_with_external_anchor.png` と同名SVGに3パネルを収録。

- **A: paired heatmap**。左に全48-feature PLIF（既存contact frequency）、右に全14-term clinical fingerprint（既存logROR）。両者の薬剤順は既存frozen Yの行順をそのまま使用：diazepam、alprazolam、triazolam、zolpidem、lorazepam、clonazepam、midazolam、temazepam、zopiclone、zaleplon。特徴量の順序も元のまま。独立の色尺度を使用し、値の再標準化・clustering・feature selectionは行わない。F01–F21はalpha1、F22–F48はalpha2で、全feature名は`data/feature_display_key.csv`とBの行ラベルに対応。NAは灰色でありゼロではない。
- **B: drug-level Spearman correlation matrix**。元runの全672相関を読み込み、再計算しない。色がrho、各セルの数字が元表のn_effective。X featureが行、clinical termが列、統計単位は薬剤で最大n=10。灰色は定数feature等でrho未定義。n_effectiveが表示されていても相関を計算可能とは限らない。これは薬剤対の45×45行列ではない。
- **C: external biological anchor**。diazepam/triazolamのNishino2008 rotarodのみ。元runのcontext-matchedな2/5 mg/kg × 15/30/60/90分の16行を使用。doseごとに別パネルとし、既存failure fractionをそのまま表示。mouse/ICR/male/oral、各報告dose/time条件n=10 mice。用量を跨ぐ平均・再正規化・n=2相関・mediation推定を行わない。元のfraction自体は前runでfailure count / sample sizeとして作成済みであり、本図作成時は再計算していない。

PLIF featureとclinical termを直接結ぶ矢印はない。drug identityはデータの共通キーであり、因果経路ではない。rotarodはexternal biological anchorであり、実証されたmediatorではない。FAERS logRORは自発報告のdisproportionalityで、incidenceやrisk ratioではない。既存Yの低count由来NA 6セルはそのまま保持。alpha2 provisional construct等の元の構造制約も残る。

## Inputs / provenance

元run: `runs/20260920_clinical_bridge_hypothesis_poc/`（commit `a39e7668ec31eae0e5fed6a0488b97d0815576f4`）。構造の元runは`runs/20260919_0724_10drug_standardized_redocking/`。

`data/`には使用したX/Y、相関表、rotarod表、clinical term mapping、freeze records、元manifest/configを**byte-identical**に保存。`input_manifest.json`に元pathとSHA256を記録。元のANALYSIS_CONTRACT.mdは`data/ANALYSIS_CONTRACT_source.md`にsnapshot。本作業は同じX↔Y比較の表示変更のみで、解析軸・入力値・feature集合を変更しない。元のrunは編集しない。

- X SHA256: `c986ba328a01fd4bc6909ca24e4724d7eea5971e5a5ea6ae0fc2faa3147b2f10`
- Y SHA256: `8475dd9e413a97075820b93531550ffa8be49dba591834dab3f4e5ccdcb01ea2`
- [Nishino2008](https://doi.org/10.1254/jphs.08107FP): 既存runで原表監査済みの数値。今回は新規抽出なし。

## Regenerate

Python、numpy、pandas、matplotlibが必要。

```bash
python runs/20260920_frozen_xy_correspondence_visualization/code/make_figure.py
```

スクリプトはfreeze hash、drug identity一致、48/14列、672セル、anchor contextを検証して描画する。相関計算ライブラリ・新規API取得・dockingを使用しない。既存data snapshotが元fileと一致しない場合は停止する。SVGはtextを保持する。PNGは180 dpi。`QC_REPORT.json`は表示用検証結果。`feature_display_key.csv`と`drug_display_order.csv`のみ表示用metadataであり新しいfingerprintではない。
