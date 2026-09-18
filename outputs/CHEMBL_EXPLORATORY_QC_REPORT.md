# Exploratory GABA-A pharmacology QC

対象activity行数: 75
対象targetは指定4組成のHomo sapiens protein complexのみ。

## Drug × receptor composition

| drug_id | alpha1/beta2/gamma2 | alpha1/beta3/gamma2 | alpha2/beta2/gamma2 | alpha2/beta3/gamma2 |
|---|---:|---:|---:|---:|
| diazepam | 12 | 15 | 1 | 16 |
| lorazepam | 0 | 0 | 0 | 0 |
| clonazepam | 0 | 0 | 0 | 0 |
| alprazolam | 1 | 0 | 1 | 0 |
| midazolam | 0 | 0 | 0 | 0 |
| temazepam | 0 | 0 | 0 | 0 |
| triazolam | 0 | 1 | 0 | 1 |
| zolpidem | 4 | 9 | 3 | 9 |
| zopiclone | 0 | 0 | 0 | 0 |
| zaleplon | 1 | 0 | 0 | 1 |

## compositionごとの薬剤数

| composition | drug count |
|---|---:|
| alpha1/beta2/gamma2 | 4 |
| alpha1/beta3/gamma2 | 3 |
| alpha2/beta2/gamma2 | 3 |
| alpha2/beta3/gamma2 | 4 |

## alpha1/alpha2共通metric-unit QC

| drug_id | common metric-unit | same beta comparison | different beta comparison |
|---|---|---:|---:|
| diazepam | Activity / % ; EC50 / nM ; Efficacy / % ; Ki / nM | 4 | 5 |
| alprazolam | Ki / nM | 1 | 0 |
| triazolam | Ki / nM | 1 | 0 |
| zolpidem | EC50 / nM ; Efficacy / % ; Ki / nM | 3 | 2 |
| zaleplon | EC50 / nM | 0 | 1 |

同じstandard_type・unitを持つ薬剤: diazepam, alprazolam, triazolam, zolpidem, zaleplon

beta subtype一致と不一致を分けた行単位QCは `results/tables/chembl_exploratory_comparability_qc.csv` に保存。standard_valueの代表値決定・平均化はしていない。

Primary datasetは `data/processed/pharmacology.csv` および `results/tables/chembl_primary_*.csv` に変更せず保持。Exploratory datasetは `data/processed/pharmacology_exploratory.csv` と `results/tables/chembl_exploratory_*.csv` に分離。
