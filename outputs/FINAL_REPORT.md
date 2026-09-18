# FINAL_REPORT: strict 3-drug GABA-A Interaction Fingerprint PoC

## Outcome

The strict pipeline was executed for diazepam, triazolam and zolpidem. ChEMBL Ki rows were retained pair-by-pair, 3D ligands were generated from the downloaded canonical SMILES, six identical-condition Vina jobs were run, and pose-level PLIP fingerprints were written where parsing succeeded.

- Pair-level Ki rows: 122; no values were averaged.
- Docking rows: 54 (up to nine poses per drug/receptor).
- PLIP rows: 238.

## Interpretation boundary

The α1 receptor is directly represented by human α1β3γ2L structure 6HUP. The α2 receptor is a provisional α2/β3/γ2 local-site construct extracted from mixed native human structure 9CTJ because an exact α2β3γ2 full pentamer with a BZD-site co-crystal ligand was not found. Accordingly, docking scores, residue fingerprints and Jaccard values are mechanistic feasibility outputs, not validated subtype comparisons.

No statistical generalization, ML, imputation, activity-type integration, unit mixing or representative-value selection was performed. Same-document/same-unit Ki ratios are kept as multiple raw pair rows.

## Pair-level Ki and IFP comparison

The following compact table shows only the highest-priority same-document Ki pairs; the complete unaveraged pair table is `strict_pharmacology_summary.csv`.

```
  drug_id alpha1_document_chembl_id  alpha1_standard_value  alpha2_standard_value standard_units  ki_alpha2_over_alpha1
 diazepam             CHEMBL1132940                   14.0                  20.00             nM               1.428571
 diazepam             CHEMBL1139460                   14.0                  20.00             nM               1.428571
 diazepam             CHEMBL1142695                   13.0                   6.60             nM               0.507692
 diazepam             CHEMBL1148356                   14.3                  26.20             nM               1.832168
 diazepam             CHEMBL1148587                   13.0                   6.60             nM               0.507692
 diazepam             CHEMBL1149126                   14.0                  20.00             nM               1.428571
 diazepam             CHEMBL1149444                   13.0                  33.00             nM               2.538462
 diazepam             CHEMBL6078658                   31.0                  22.00             nM               0.709677
triazolam             CHEMBL1132940                    0.8                   0.59             nM               0.737500
 zolpidem             CHEMBL1131013                   26.7                 156.00             nM               5.842697
 zolpidem             CHEMBL1132940                   26.7                 156.00             nM               5.842697
 zolpidem             CHEMBL1139460                   27.0                 160.00             nM               5.925926
 zolpidem             CHEMBL1142695                   27.0                 103.00             nM               3.814815
 zolpidem             CHEMBL1148356                   52.1                 255.10             nM               4.896353
 zolpidem             CHEMBL1149126                   27.0                 160.00             nM               5.925926
 zolpidem             CHEMBL1268934                   26.7                 156.00             nM               5.842697
```

```
  drug_id  jaccard_alpha1_vs_alpha2  changed_interaction_feature_count  alpha1_features  alpha2_features fingerprint_status
 diazepam                  0.294118                                 12               13                9 computed_from_PLIP
triazolam                  0.111111                                 16               11                9 computed_from_PLIP
 zolpidem                  0.285714                                 10               11                7 computed_from_PLIP
```

## Files

- `strict_pharmacology_summary.csv`: raw pair-level Ki comparisons and α2/α1 ratios.
- `docking_results.csv`: drug/receptor/pose IDs and Vina scores.
- `residue_mapping.csv`: mapping audit including unresolved positions.
- `fingerprint_binary.csv`, `fingerprint_frequency.csv`: pose-level IFP outputs.
- `alpha1_alpha2_similarity.csv`: Jaccard and changed-feature counts.
- `QC_REPORT.md`: explicit structural, box, PDBQT and interaction QC.

---

---

---

---

---

---

---

---

# FOUR_DRUG_POC_REPORT

This is a four-drug Tier-1 mechanistic feasibility extension. Existing strict three-drug results remain in their original files.

## Primary pharmacology

Tier-1 raw pair rows: **24**; summary groups: **24**. Ki-only exact comparison-context rows retained for integration: **18** (no cross-reference representative was selected).

   drug_id  alpha1_composition  alpha2_composition metric unit alpha1_reference_id alpha2_reference_id  alpha1_assay  alpha2_assay  number_of_comparable_pairs  alpha1_value_median  alpha2_value_median  alpha2_over_alpha1_ratio_median  log10_ratio_median     reference      assay_id  reference_count  assay_count                                   comparison_context_id                                                                                  selection_rule  jaccard_similarity  shared_features  alpha1_only_features  alpha2_only_features  changed_features  cosine_similarity_frequency  euclidean_distance_frequency  valid_features  unresolved_mapping_included
alprazolam alpha1/beta2/gamma2 alpha2/beta2/gamma2     Ki   nM       CHEMBL5143601       CHEMBL5143601 CHEMBL5146040 CHEMBL5146041                           1                  0.8                 0.60                         0.750000           -0.124939 CHEMBL5143601 CHEMBL5146040                1            1 CHEMBL5143601|CHEMBL5143601|CHEMBL5146040|CHEMBL5146041 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.777778                7                     0                     2                 2                     0.812721                      0.666667               9                         True
  diazepam alpha1/beta2/gamma2 alpha2/beta2/gamma2     Ki   nM       CHEMBL5143601       CHEMBL5143601 CHEMBL5146040 CHEMBL5146041                           1                 14.0                20.00                         1.428571            0.154902 CHEMBL5143601 CHEMBL5146040                1            1 CHEMBL5143601|CHEMBL5143601|CHEMBL5146040|CHEMBL5146041 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1132940       CHEMBL1132940  CHEMBL678329  CHEMBL681248                           1                 14.0                20.00                         1.428571            0.154902 CHEMBL1132940  CHEMBL678329                1            1   CHEMBL1132940|CHEMBL1132940|CHEMBL678329|CHEMBL681248 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1139460       CHEMBL1139460  CHEMBL919087  CHEMBL919088                           1                 14.0                20.00                         1.428571            0.154902 CHEMBL1139460  CHEMBL919087                1            1   CHEMBL1139460|CHEMBL1139460|CHEMBL919087|CHEMBL919088 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1142695       CHEMBL1142695  CHEMBL829006  CHEMBL829005                           1                 13.0                 6.60                         0.507692           -0.294399 CHEMBL1142695  CHEMBL829006                1            1   CHEMBL1142695|CHEMBL1142695|CHEMBL829006|CHEMBL829005 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1148356       CHEMBL1148356  CHEMBL677034  CHEMBL682709                           1                 14.3                26.20                         1.832168            0.262965 CHEMBL1148356  CHEMBL677034                1            1   CHEMBL1148356|CHEMBL1148356|CHEMBL677034|CHEMBL682709 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1148587       CHEMBL1148587  CHEMBL864002  CHEMBL864003                           1                 13.0                 6.60                         0.507692           -0.294399 CHEMBL1148587  CHEMBL864002                1            1   CHEMBL1148587|CHEMBL1148587|CHEMBL864002|CHEMBL864003 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1149126       CHEMBL1149126  CHEMBL676826  CHEMBL681841                           1                 14.0                20.00                         1.428571            0.154902 CHEMBL1149126  CHEMBL676826                1            1   CHEMBL1149126|CHEMBL1149126|CHEMBL676826|CHEMBL681841 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1149444       CHEMBL1149444  CHEMBL853068  CHEMBL853070                           1                 13.0                33.00                         2.538462            0.404571 CHEMBL1149444  CHEMBL853068                1            1   CHEMBL1149444|CHEMBL1149444|CHEMBL853068|CHEMBL853070 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
  diazepam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL6078658       CHEMBL6078658 CHEMBL6080337 CHEMBL6080338                           1                 31.0                22.00                         0.709677           -0.148939 CHEMBL6078658 CHEMBL6080337                1            1 CHEMBL6078658|CHEMBL6078658|CHEMBL6080337|CHEMBL6080338 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
 triazolam alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1132940       CHEMBL1132940  CHEMBL678329  CHEMBL681248                           1                  0.8                 0.59                         0.737500           -0.132238 CHEMBL1132940  CHEMBL678329                1            1   CHEMBL1132940|CHEMBL1132940|CHEMBL678329|CHEMBL681248 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.222222                2                     0                     7                 7                     0.503739                      0.955814               9                         True
  zolpidem alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1131013       CHEMBL1131013  CHEMBL822162  CHEMBL824294                           1                 26.7               156.00                         5.842697            0.766613 CHEMBL1131013  CHEMBL822162                1            1   CHEMBL1131013|CHEMBL1131013|CHEMBL822162|CHEMBL824294 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True
  zolpidem alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1132940       CHEMBL1132940  CHEMBL678329  CHEMBL681248                           1                 26.7               156.00                         5.842697            0.766613 CHEMBL1132940  CHEMBL678329                1            1   CHEMBL1132940|CHEMBL1132940|CHEMBL678329|CHEMBL681248 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True
  zolpidem alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1139460       CHEMBL1139460  CHEMBL919087  CHEMBL919088                           1                 27.0               160.00                         5.925926            0.772756 CHEMBL1139460  CHEMBL919087                1            1   CHEMBL1139460|CHEMBL1139460|CHEMBL919087|CHEMBL919088 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True
  zolpidem alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1142695       CHEMBL1142695  CHEMBL829006  CHEMBL829005                           1                 27.0               103.00                         3.814815            0.581473 CHEMBL1142695  CHEMBL829006                1            1   CHEMBL1142695|CHEMBL1142695|CHEMBL829006|CHEMBL829005 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True
  zolpidem alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1148356       CHEMBL1148356  CHEMBL677034  CHEMBL682709                           1                 52.1               255.10                         4.896353            0.689873 CHEMBL1148356  CHEMBL677034                1            1   CHEMBL1148356|CHEMBL1148356|CHEMBL677034|CHEMBL682709 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True
  zolpidem alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1149126       CHEMBL1149126  CHEMBL676826  CHEMBL681841                           1                 27.0               160.00                         5.925926            0.772756 CHEMBL1149126  CHEMBL676826                1            1   CHEMBL1149126|CHEMBL1149126|CHEMBL676826|CHEMBL681841 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True
  zolpidem alpha1/beta3/gamma2 alpha2/beta3/gamma2     Ki   nM       CHEMBL1268934       CHEMBL1268934 CHEMBL1273615 CHEMBL1273616                           1                 26.7               156.00                         5.842697            0.766613 CHEMBL1268934 CHEMBL1273615                1            1 CHEMBL1268934|CHEMBL1268934|CHEMBL1273615|CHEMBL1273616 Ki-only; exact composition/metric/unit/reference/assay context; median only within this context            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True

## IFP results

   drug_id  jaccard_similarity  shared_features  alpha1_only_features  alpha2_only_features  changed_features  cosine_similarity_frequency  euclidean_distance_frequency  valid_features  unresolved_mapping_included
  diazepam            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True
alprazolam            0.777778                7                     0                     2                 2                     0.812721                      0.666667               9                         True
 triazolam            0.222222                2                     0                     7                 7                     0.503739                      0.955814               9                         True
  zolpidem            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True

Alprazolam docking rows: **18**; PLIP interaction rows in four-drug set: **339**; unresolved mapping rows: **60**.

Differential mapped-feature counts: {'alprazolam': 9, 'diazepam': 8, 'triazolam': 9, 'zolpidem': 6}.

## Drug-specific differential features

                   feature changed_in_drugs  n_changed_drugs  drug_specific drug_specific_to
 BZD_SITE_193|halogen_bond        triazolam                1           True        triazolam
BZD_SITE_205|hydrogen_bond         diazepam                1           True         diazepam

## Comparison with the previous three-drug PoC

  drug_id  jaccard_similarity  shared_features  alpha1_only_features  alpha2_only_features  changed_features  cosine_similarity_frequency  euclidean_distance_frequency  valid_features  unresolved_mapping_included  jaccard_alpha1_vs_alpha2  changed_interaction_feature_count  alpha1_features  alpha2_features fingerprint_status
 diazepam            0.555556                5                     0                     4                 4                     0.787044                      0.555556               9                         True                  0.294118                                 12               13                9 computed_from_PLIP
triazolam            0.222222                2                     0                     7                 7                     0.503739                      0.955814               9                         True                  0.111111                                 16               11                9 computed_from_PLIP
 zolpidem            0.571429                4                     0                     3                 3                     0.862796                      0.496904               7                         True                  0.285714                                 10               11                7 computed_from_PLIP

Alprazolam expands the feature space and provides a fourth qualitative comparison. Any apparent correspondence between Ki ratios and IFP similarity is descriptive only; n=4 is insufficient for a generalization or significance claim.

## Answers to the requested interpretation questions

1. **Pharmacology difference:** alprazolam: [0.75 (log10=-0.125; Ki/nM; 1 pair)] across 1 exact context(s); diazepam: [1.43 (log10=0.155; Ki/nM; 1 pair), 1.43 (log10=0.155; Ki/nM; 1 pair), 1.43 (log10=0.155; Ki/nM; 1 pair), 0.508 (log10=-0.294; Ki/nM; 1 pair), 1.83 (log10=0.263; Ki/nM; 1 pair), 0.508 (log10=-0.294; Ki/nM; 1 pair), 1.43 (log10=0.155; Ki/nM; 1 pair), 2.54 (log10=0.405; Ki/nM; 1 pair), 0.71 (log10=-0.149; Ki/nM; 1 pair)] across 9 exact context(s); triazolam: [0.737 (log10=-0.132; Ki/nM; 1 pair)] across 1 exact context(s); zolpidem: [5.84 (log10=0.767; Ki/nM; 1 pair), 5.84 (log10=0.767; Ki/nM; 1 pair), 5.93 (log10=0.773; Ki/nM; 1 pair), 3.81 (log10=0.581; Ki/nM; 1 pair), 4.9 (log10=0.69; Ki/nM; 1 pair), 5.93 (log10=0.773; Ki/nM; 1 pair), 5.84 (log10=0.767; Ki/nM; 1 pair)] across 7 exact context(s). Ratios are listed per exact reference/assay context; no values from different references were averaged. Alprazolam has only alpha1/beta2/gamma2 versus alpha2/beta2/gamma2 Tier-1 evidence; diazepam also has one beta2 context alongside its beta3 contexts. These composition differences are flagged for review against the fixed beta3 docking background.
2. **IFP similarity:** diazepam: Jaccard=0.556, cosine=0.787, changed=4; alprazolam: Jaccard=0.778, cosine=0.813, changed=2; triazolam: Jaccard=0.222, cosine=0.504, changed=7; zolpidem: Jaccard=0.571, cosine=0.863, changed=3. Similarity uses only features with resolved mapping on both receptor sides; unresolved rows are excluded from the metric and shown as gray in the heatmaps.
3. **Apparent pharmacology–IFP trend:** the four ratios and similarities are not monotonic (for example zolpidem has the largest Ki ratio but not the lowest Jaccard). No general trend is claimed from n=4.
4. **Drug-specific differences:** resolved drug-specific features are BZD_SITE_193|halogen_bond (triazolam), BZD_SITE_205|hydrogen_bond (diazepam). They are descriptive candidates, not biomarkers.
5. **Zolpidem:** resolved differential features are BZD_SITE_100|hydrophobic_interaction [alpha2_only], BZD_SITE_102|hydrophobic_interaction [alpha2_frequency_higher], BZD_SITE_103|hydrophobic_interaction [alpha1_frequency_higher], BZD_SITE_193|hydrophobic_interaction [alpha2_only], BZD_SITE_203|hydrophobic_interaction [alpha1_frequency_higher], BZD_SITE_210|hydrophobic_interaction [alpha2_only]. The drug-specific table does not identify a feature unique to zolpidem; unresolved mapping features are not interpreted.
6. **Effect of adding alprazolam:** the common feature space now contains four drugs and the four-drug similarity table adds alprazolam (Jaccard 0.778; changed features 2). The prior three-drug rows and files remain preserved, while the beta2 pharmacology versus beta3 docking mismatch limits direct mechanistic integration for alprazolam.
7. **Largest confounder:** the provisional alpha2 local construct/chain background and transferred box are not an exact matched alpha1/alpha2 structural pair; pharmacology contexts that use beta2 while docking uses beta3 (all alprazolam contexts and one diazepam context) are additional composition confounders.
8. **Value of more drugs:** additional compounds are useful for testing whether resolved feature-level patterns recur, provided they have same-metric, same-unit, composition-matched Tier-1 pharmacology and comparable structures. This PoC does not estimate predictive performance.

## QC and limitations

The existing 6HUP alpha1 receptor, provisional 9CTJ-derived alpha2 local construct, transferred DZP box, receptor preparation, residue mapping and Vina parameters were reused. The alpha2 construct is not an exact alpha2beta3gamma2 full pentamer, and the two receptor backgrounds are not fully identical. Gray heatmap cells indicate unresolved mapping; zero indicates no interaction observed only when pose data and mapping were available. GtoPdb was not integrated in this phase.

The archived ChEMBL table does not retain assay_description for all Tier-1 rows, so annotation-versus-assay-description contradiction checking is reported as unavailable rather than inferred.

The exact strict three-drug files under `outputs/` were left intact. Four-drug combined tables use `_4drug` names; the required four-drug `data/processed/fingerprint_binary.csv` and `fingerprint_frequency.csv` are generated separately from those strict outputs.

Four-drug combined tables: `results/tables/docking_results_4drug.csv` (72 pose rows; 9 poses per drug/receptor), `results/tables/plip_interactions_4drug.csv`, `results/tables/fingerprint_binary_4drug.csv`, and `results/tables/fingerprint_frequency_4drug.csv`. A computational/biological QC checklist is in `outputs/FOUR_DRUG_QC_REPORT.md`.

