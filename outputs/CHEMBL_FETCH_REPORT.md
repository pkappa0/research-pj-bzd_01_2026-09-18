# ChEMBL pharmacology取得結果

- 対象薬剤: 10
- exact preferred-name match: 10
- exact match失敗: 0
- 取得activity行数: 434
- target候補数: 90
- API metadata: `data/raw/pharmacology/chembl/chembl_api_metadata.json`

## 分子情報

| input drug name | ChEMBL ID | canonical SMILES | InChIKey |
|---|---|---|---|
| Diazepam | CHEMBL12 | `CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21` | AAOVKJBEBIDNHE-UHFFFAOYSA-N |
| Lorazepam | CHEMBL580 | `O=C1Nc2ccc(Cl)cc2C(c2ccccc2Cl)=NC1O` | DIWRORZWFLOCLC-UHFFFAOYSA-N |
| Clonazepam | CHEMBL452 | `O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1` | DGBIGWXXNGSACT-UHFFFAOYSA-N |
| Alprazolam | CHEMBL661 | `Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2` | VREFGVBLTWBCJP-UHFFFAOYSA-N |
| Midazolam | CHEMBL655 | `Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2` | DDLIGBOFAVUZHB-UHFFFAOYSA-N |
| Temazepam | CHEMBL967 | `CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21` | SEQDDYPDSLOBDC-UHFFFAOYSA-N |
| Triazolam | CHEMBL646 | `Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2` | JOFWLTCLBGQGBO-UHFFFAOYSA-N |
| Zolpidem | CHEMBL911 | `Cc1ccc(-c2nc3ccc(C)cn3c2CC(=O)N(C)C)cc1` | ZAFYATHCZYHLPB-UHFFFAOYSA-N |
| Zopiclone | CHEMBL135400 | `CN1CCN(C(=O)OC2c3nccnc3C(=O)N2c2ccc(Cl)cn2)CC1` | GBBSUAFBMRNDJC-UHFFFAOYSA-N |
| Zaleplon | CHEMBL1521 | `CCN(C(C)=O)c1cccc(-c2ccnc3c(C#N)cnn23)c1` | HUNXMJYCHXQEGX-UHFFFAOYSA-N |

## 薬剤別activity件数

| drug_id | input drug name | activity count |
|---|---|---|
| diazepam | Diazepam | 263 |
| lorazepam | Lorazepam | 2 |
| clonazepam | Clonazepam | 16 |
| alprazolam | Alprazolam | 16 |
| midazolam | Midazolam | 3 |
| temazepam | Temazepam | 3 |
| triazolam | Triazolam | 12 |
| zolpidem | Zolpidem | 109 |
| zopiclone | Zopiclone | 4 |
| zaleplon | Zaleplon | 6 |

## target別件数

| target ChEMBL ID | target name | target class | receptor subtype | count |
|---|---|---|---|---|
| CHEMBL1907607 | GABA-A receptor; anion channel | generic_gabaa_target | nan | 104 |
| CHEMBL2094130 | GABA-A receptor; alpha-2/beta-3/gamma-2 | composite_target | alpha2/beta3/gamma2 | 27 |
| CHEMBL2094120 | GABA-A receptor; alpha-3/beta-3/gamma-2 | composite_target | alpha3/beta3/gamma2 | 26 |
| CHEMBL2094121 | GABA-A receptor; alpha-1/beta-3/gamma-2 | composite_target | alpha1/beta3/gamma2 | 25 |
| CHEMBL343 | Gamma-aminobutyric acid receptor subunit alpha-1 | subunit_specific | alpha1 | 22 |
| CHEMBL2094107 | GABA-A receptor; anion channel | generic_gabaa_target | nan | 21 |
| CHEMBL2095167 | GABA-A receptor; alpha-1/beta-2/gamma-2 | composite_target | alpha1/beta2/gamma2 | 21 |
| CHEMBL2094122 | GABA-A receptor; alpha-5/beta-3/gamma-2 | composite_target | alpha5/beta3/gamma2 | 20 |
| CHEMBL2095172 | GABA-A receptor; alpha-1/beta-2/gamma-2 | composite_target | alpha1/beta2/gamma2 | 18 |
| CHEMBL2093872 | GABA-A receptor; anion channel | generic_gabaa_target | nan | 16 |
| CHEMBL1962 | Gamma-aminobutyric acid receptor subunit alpha-1 | subunit_specific | alpha1 | 14 |
| CHEMBL2111365 | GABA A receptor alpha-6/beta-2/gamma-2 | composite_target | alpha6/beta2/gamma2 | 11 |
| CHEMBL2111343 | GABA A receptor alpha-3/beta-2/gamma-2 | composite_target | alpha3/beta2/gamma2 | 9 |
| CHEMBL2111374 | GABA A receptor alpha-5/beta-3/gamma-2 | composite_target | alpha5/beta3/gamma2 | 9 |
| CHEMBL5112 | Gamma-aminobutyric acid receptor subunit alpha-5 | subunit_specific | alpha5 | 8 |
| CHEMBL2111339 | GABA A receptor alpha-3/beta-2/gamma-2 | composite_target | alpha3/beta2/gamma2 | 7 |
| CHEMBL3026 | Gamma-aminobutyric acid receptor subunit alpha-3 | subunit_specific | alpha3 | 7 |
| CHEMBL341 | Gamma-aminobutyric acid receptor subunit alpha-2 | subunit_specific | alpha2 | 7 |
| CHEMBL4296062 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-3/gamma-2 | composite_target | alpha1/beta3/gamma2 | 7 |
| CHEMBL4956 | Gamma-aminobutyric acid receptor subunit alpha-2 | subunit_specific | alpha2 | 6 |
| CHEMBL2094133 | GABA-A receptor; anion channel | generic_gabaa_target | nan | 5 |
| CHEMBL2095190 | GABA-A receptor; alpha-6/beta-3/gamma-2 | composite_target | alpha6/beta3/gamma2 | 5 |
| CHEMBL2111327 | GABA A receptor alpha-2/beta-2/gamma-2 | composite_target | alpha2/beta2/gamma2 | 5 |
| CHEMBL2111413 | GABA A receptor alpha-2/beta-2/gamma-2 | composite_target | alpha2/beta2/gamma2 | 5 |
| CHEMBL5291947 | Gamma-aminobutyric acid receptor subunit alpha-1/alpha-6/beta-3/gamma-2 | composite_target | alpha1/alpha6/beta3/gamma2 | 5 |
| CHEMBL4296054 | Gamma-aminobutyric acid receptor subunit alpha-6/beta-3/gamma-2 | composite_target | alpha6/beta3/gamma2 | 4 |
| CHEMBL2111370 | GABA A receptor alpha-6/beta-2/gamma-2 | composite_target | alpha6/beta2/gamma2 | 3 |
| CHEMBL2111392 | GABA A receptor alpha-1/beta-1/gamma-2 | composite_target | alpha1/beta1/gamma2 | 2 |
| CHEMBL300 | Gamma-aminobutyric acid receptor subunit alpha-5 | subunit_specific | alpha5 | 2 |
| CHEMBL328 | Gamma-aminobutyric acid receptor subunit alpha-3 | subunit_specific | alpha3 | 2 |
| CHEMBL4296052 | Gamma-aminobutyric acid receptor subunit alpha-6/beta-3 | composite_target | alpha6/beta3 | 2 |
| CHEMBL4523641 | Gamma-aminobutyric acid receptor subunit alpha-5/beta-2/gamma-2 | composite_target | alpha5/beta2/gamma2 | 2 |
| CHEMBL2472 | Gamma-aminobutyric acid receptor subunit alpha-4 | subunit_specific | alpha4 | 1 |
| CHEMBL296 | Gamma-aminobutyric acid receptor subunit gamma-1 | subunit_specific | gamma1 | 1 |
| CHEMBL3883322 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-2/beta-3/gamma-2 | composite_target | alpha1/beta2/beta3/gamma2 | 1 |
| CHEMBL3885577 | Gamma-aminobutyric acid receptor subunit alpha-5/beta-3/gamma-3 | composite_target | alpha5/beta3/gamma3 | 1 |
| CHEMBL4296047 | Gamma-aminobutyric acid receptor subunit alpha-3/beta-3/gamma-2 | composite_target | alpha3/beta3/gamma2 | 1 |
| CHEMBL5291948 | Gamma-aminobutyric acid receptor subunit alpha-6/beta-3/alpha-1 | composite_target | alpha6/beta3/alpha1 | 1 |
| CHEMBL5291950 | Gamma-aminobutyric acid receptor subunit gamma-2/beta-3 | composite_target | gamma2/beta3 | 1 |

## standard_type別件数

| standard_type | count |
|---|---|
| Ki | 188 |
| IC50 | 63 |
| Activity | 48 |
| EC50 | 26 |
| Efficacy | 25 |
| AC50 | 20 |
| Ratio | 11 |
| Inhibition | 10 |
| GABA ratio | 8 |
| GABA shift | 5 |
| Kd | 4 |
| Displacement | 3 |
| ED50 | 2 |
| ID50 | 2 |
| TBPS shift | 2 |
| Change in Cl- current | 2 |
| Cl - current change | 2 |
| GS | 2 |
| Control | 2 |
| Clonazepam | 1 |
| FC | 1 |
| GR | 1 |
| Log IC50 | 1 |
| Binding affinity | 1 |
| Shift | 1 |
| GABA current | 1 |
| [35S] TBPS shift | 1 |
| max activation | 1 |

## target分類

| class | count |
|---|---|
| composite_target | 218 |
| generic_gabaa_target | 146 |
| subunit_specific | 70 |

## α1/α2候補

全候補を保存し、採用可否は決定していません。activityが存在する候補には件数を付けています。

| target ID | target name | type | class | subtype | activity count |
|---|---|---|---|---|---|
| CHEMBL1962 | Gamma-aminobutyric acid receptor subunit alpha-1 | SINGLE PROTEIN | subunit_specific | alpha1 | 14 |
| CHEMBL3139 | Gamma-aminobutyric acid receptor subunit alpha-1 | SINGLE PROTEIN | subunit_specific | alpha1 | 0 |
| CHEMBL4956 | Gamma-aminobutyric acid receptor subunit alpha-2 | SINGLE PROTEIN | subunit_specific | alpha2 | 6 |
| CHEMBL341 | Gamma-aminobutyric acid receptor subunit alpha-2 | SINGLE PROTEIN | subunit_specific | alpha2 | 7 |
| CHEMBL343 | Gamma-aminobutyric acid receptor subunit alpha-1 | SINGLE PROTEIN | subunit_specific | alpha1 | 22 |
| CHEMBL3883322 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-2/beta-3/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta2/beta3/gamma2 | 1 |
| CHEMBL3885570 | Gamma-aminobutyric acid receptor subunit alpha-1/ beta-1 | PROTEIN COMPLEX | composite_target | alpha1/beta1 | 0 |
| CHEMBL3885571 | Gamma-aminobutyric acid receptor subunit alpha-2/beta-2 | PROTEIN COMPLEX | composite_target | alpha2/beta2 | 0 |
| CHEMBL4296044 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-2/delta | PROTEIN COMPLEX | composite_target | alpha1/beta2/delta | 0 |
| CHEMBL4296048 | Gamma-aminobutyric acid receptor subunit alpha-2/beta-3/gamma-2 | PROTEIN COMPLEX | composite_target | alpha2/beta3/gamma2 | 0 |
| CHEMBL4296058 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-3/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta3/gamma2 | 0 |
| CHEMBL4296059 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-2/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta2/gamma2 | 0 |
| CHEMBL4296060 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-3 | PROTEIN COMPLEX | composite_target | alpha1/beta3 | 0 |
| CHEMBL4296061 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-3/delta | PROTEIN COMPLEX | composite_target | alpha1/beta3/delta | 0 |
| CHEMBL4296062 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-3/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta3/gamma2 | 7 |
| CHEMBL4296063 | Gamma-aminobutyric acid receptor subunit alpha-1/beta-2 | PROTEIN COMPLEX | composite_target | alpha1/beta2 | 0 |
| CHEMBL4296064 | Gamma-aminobutyric acid receptor subunit alpha-2/beta-3 | PROTEIN COMPLEX | composite_target | alpha2/beta3 | 0 |
| CHEMBL4523638 | Gamma-aminobutyric acid receptor subunit alpha-1/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/gamma2 | 0 |
| CHEMBL5291947 | Gamma-aminobutyric acid receptor subunit alpha-1/alpha-6/beta-3/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/alpha6/beta3/gamma2 | 5 |
| CHEMBL5291948 | Gamma-aminobutyric acid receptor subunit alpha-6/beta-3/alpha-1 | PROTEIN COMPLEX | composite_target | alpha6/beta3/alpha1 | 1 |
| CHEMBL5291949 | Gamma-aminobutyric acid receptor subunit alpha-2/beta-1/gamma-2 | PROTEIN COMPLEX | composite_target | alpha2/beta1/gamma2 | 0 |
| CHEMBL5303741 | Gamma-aminobutyric acid receptor subunit alpha-1/alpha-2/beta-2/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/alpha2/beta2/gamma2 | 0 |
| CHEMBL1907597 | GABA-A receptor; GABA-A site (alpha1/beta2 interface) | PROTEIN COMPLEX | composite_target | alpha1/beta2 | 0 |
| CHEMBL2094121 | GABA-A receptor; alpha-1/beta-3/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta3/gamma2 | 25 |
| CHEMBL2094130 | GABA-A receptor; alpha-2/beta-3/gamma-2 | PROTEIN COMPLEX | composite_target | alpha2/beta3/gamma2 | 27 |
| CHEMBL2095167 | GABA-A receptor; alpha-1/beta-2/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta2/gamma2 | 21 |
| CHEMBL2095172 | GABA-A receptor; alpha-1/beta-2/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta2/gamma2 | 18 |
| CHEMBL4106151 | GABA-A receptor alpha-1/beta-3 | PROTEIN COMPLEX | composite_target | alpha1/beta3 | 0 |
| CHEMBL4680049 | GABA-A receptor alpha-1/beta-1 | PROTEIN COMPLEX | composite_target | alpha1/beta1 | 0 |
| CHEMBL2111413 | GABA A receptor alpha-2/beta-2/gamma-2 | PROTEIN COMPLEX | composite_target | alpha2/beta2/gamma2 | 5 |
| CHEMBL2111327 | GABA A receptor alpha-2/beta-2/gamma-2 | PROTEIN COMPLEX | composite_target | alpha2/beta2/gamma2 | 5 |
| CHEMBL2111392 | GABA A receptor alpha-1/beta-1/gamma-2 | PROTEIN COMPLEX | composite_target | alpha1/beta1/gamma2 | 2 |

## 注意

- Ki/Kd/IC50/EC50/AC50等は統合・平均化していません。
- `generic_gabaa_target`、`composite_target`、`subunit_specific`を別分類で保持しています。
- raw JSONL/API pageは `data/raw/pharmacology/chembl/` に保存しています。
- exact match失敗薬剤はありませんでした。
- 一部standard_typeはKi/Kd/IC50/EC50以外であり、後段で採用基準を決定してください。
