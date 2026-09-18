# Primary α1β3γ2 vs α2β3γ2 QC

対象target:

- `CHEMBL2094121`: GABA-A receptor; alpha-1/beta-3/gamma-2 (Homo sapiens)
- `CHEMBL2094130`: GABA-A receptor; alpha-2/beta-3/gamma-2 (Homo sapiens)

raw primary activity rows: 52
drug × target × standard_typeの集計は別CSVに保存

## Drug × target

| drug_id | α1β3γ2 | α2β3γ2 |
|---|---:|---:|
| diazepam | 15 | 16 |
| lorazepam | 0 | 0 |
| clonazepam | 0 | 0 |
| alprazolam | 0 | 0 |
| midazolam | 0 | 0 |
| temazepam | 0 | 0 |
| triazolam | 1 | 1 |
| zolpidem | 9 | 9 |
| zopiclone | 0 | 0 |
| zaleplon | 0 | 1 |

## 両targetにactivityがある薬剤

- `diazepam`
- `triazolam`
- `zolpidem`

## α1のみ / α2のみ / どちらも無い薬剤

α1のみ: なし
α2のみ: zaleplon
どちらも無し: lorazepam, clonazepam, alprazolam, midazolam, temazepam, zopiclone

## 同じ指標・同じ単位で比較可能な薬剤

3剤。ここでの比較可能性は、両targetに同じ `standard_type` と同じ非欠損 `standard_units` の行が少なくとも1つあることだけで判定しています。値の選択、平均、assay間統合はしていません。

## 注意

- generic GABA-A target、single-protein α1/α2 targetはこのprimary表に含めていません。
- Ki/Kd/IC50/EC50等、assay、document、relation、unit、pChEMBLは行単位で保持しています。
- 複数値の分布は `chembl_primary_duplicate_value_distribution.csv` に保存し、代表値は決定していません。
- 全表は `results/tables/chembl_primary_*.csv` に保存しています。
