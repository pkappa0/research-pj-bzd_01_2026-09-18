Run ID: 20260917_1604_multireceptor_feature
Parent Run: legacy_unversioned_4drug_poc
Analysis Phase: multi-receptor fingerprint feature engineering
Started: 2026-09-17T16:04:41+09:00
Completed: 2026-09-17T16:04:57+09:00

# MULTIRECEPTOR_FINGERPRINT_REPORT

This phase is a multi-receptor interaction fingerprint feature-engineering PoC. Primary analysis is drug-to-drug comparison within a fixed receptor; alpha1-vs-alpha2 similarity remains secondary/QC.

## Dataset and primary feature construction

Four drugs: diazepam, alprazolam, triazolam, zolpidem. Alpha1 frequency feature count: **11**; alpha2 frequency feature count: **9**; combined feature columns: **200**.
PLIP interaction rows: **339**; unresolved interaction rows excluded from primary: **177**; observed interaction types: halogen_bond, hydrogen_bond, hydrophobic_interaction, pi_stack; observed distance fields: {'dist': 285, 'centdist': 35, 'dist_d-a': 19}.
Distance extraction uses PLIP dist for hydrophobic/halogen, dist-d-a for hydrogen bonds, and centdist for pi stacking. n_observed is the number of poses with the interaction; distance_n_observed counts numeric distance observations. No distance value was imputed.

## Primary within-receptor drug similarity

### Alpha1

receptor     drug_a     drug_b  jaccard_similarity_binary  shared_features  union_features  cosine_similarity_frequency  euclidean_distance_frequency  n_frequency_features
  alpha1   diazepam alprazolam                   0.636364                7              11                     0.908922                      0.400617                    11
  alpha1   diazepam  triazolam                   0.333333                3               9                     0.747666                      0.666667                    11
  alpha1   diazepam   zolpidem                   0.272727                3              11                     0.696607                      0.745356                    11
  alpha1 alprazolam  triazolam                   0.400000                4              10                     0.765118                      0.657342                    11
  alpha1 alprazolam   zolpidem                   0.600000                6              10                     0.740822                      0.702728                    11
  alpha1  triazolam   zolpidem                   0.428571                3               7                     0.877876                      0.484322                    11

### Alpha2

receptor     drug_a     drug_b  jaccard_similarity_binary  shared_features  union_features  cosine_similarity_frequency  euclidean_distance_frequency  n_frequency_features
  alpha2   diazepam alprazolam                   0.777778                7               9                     0.947331                      0.458123                     9
  alpha2   diazepam  triazolam                   0.444444                4               9                     0.652192                      0.753592                     9
  alpha2   diazepam   zolpidem                   0.555556                5               9                     0.874028                      0.415740                     9
  alpha2 alprazolam  triazolam                   0.625000                5               8                     0.772410                      0.728604                     9
  alpha2 alprazolam   zolpidem                   0.750000                6               8                     0.871365                      0.577350                     9
  alpha2  triazolam   zolpidem                   0.571429                4               7                     0.753094                      0.647884                     9

Pairwise summary:
alpha1: diazepam-alprazolam: J=0.636, cosine=0.909, Euclidean=0.401; diazepam-triazolam: J=0.333, cosine=0.748, Euclidean=0.667; diazepam-zolpidem: J=0.273, cosine=0.697, Euclidean=0.745; alprazolam-triazolam: J=0.400, cosine=0.765, Euclidean=0.657; alprazolam-zolpidem: J=0.600, cosine=0.741, Euclidean=0.703; triazolam-zolpidem: J=0.429, cosine=0.878, Euclidean=0.484
alpha2: diazepam-alprazolam: J=0.778, cosine=0.947, Euclidean=0.458; diazepam-triazolam: J=0.444, cosine=0.652, Euclidean=0.754; diazepam-zolpidem: J=0.556, cosine=0.874, Euclidean=0.416; alprazolam-triazolam: J=0.625, cosine=0.772, Euclidean=0.729; alprazolam-zolpidem: J=0.750, cosine=0.871, Euclidean=0.577; triazolam-zolpidem: J=0.571, cosine=0.753, Euclidean=0.648

## Drug-discriminating features

### Alpha1 highest-variance rows

                                    feature     statistic  variance  drug_high   drug_low
          alpha1_BZD_SITE_156_hydrogen_bond mean_distance  0.098911   zolpidem alprazolam
alpha1_BZD_SITE_102_hydrophobic_interaction mean_distance  0.061733   zolpidem alprazolam
alpha1_BZD_SITE_100_hydrophobic_interaction     frequency  0.052469   diazepam  triazolam
alpha1_BZD_SITE_100_hydrophobic_interaction mean_distance  0.043759 alprazolam   diazepam
alpha1_BZD_SITE_156_hydrophobic_interaction mean_distance  0.043512   zolpidem alprazolam

### Alpha2 highest-variance rows

                                    feature     statistic  variance  drug_high  drug_low
               alpha2_BZD_SITE_102_pi_stack mean_distance  0.100128 alprazolam  diazepam
alpha2_BZD_SITE_102_hydrophobic_interaction     frequency  0.057613  triazolam  diazepam
alpha2_BZD_SITE_156_hydrophobic_interaction     frequency  0.037037 alprazolam  diazepam
               alpha2_BZD_SITE_102_pi_stack     frequency  0.027778 alprazolam triazolam
alpha2_BZD_SITE_210_hydrophobic_interaction mean_distance  0.022008   zolpidem  diazepam

Variance/range rankings are descriptive only (n=4).

## Pharmacology label candidates

   drug_id  alpha1_Ki  alpha2_Ki  alpha2_over_alpha1_ratio  log10_alpha2_over_alpha1_ratio alpha1_Ki_raw_context_values   alpha2_Ki_raw_context_values                                                   ratio_raw_context_values  comparison_context_count                                                                  composition_contexts                                                                        aggregation_rule                       label_status
  diazepam       14.0      20.00                  1.428571                        0.154902 14;14;14.3;13;13;13;14;14;31 20;20;26.2;6.6;6.6;33;20;20;22 1.42857;1.42857;1.83217;0.507692;0.507692;2.53846;1.42857;1.42857;0.709677                         9 alpha1/beta2/gamma2 vs alpha2/beta2/gamma2;alpha1/beta3/gamma2 vs alpha2/beta3/gamma2 median across exact Tier-1 Ki contexts for label candidate; raw context values retained candidate_only_not_primary_outcome
alprazolam        0.8       0.60                  0.750000                       -0.124939                          0.8                            0.6                                                                       0.75                         1                                            alpha1/beta2/gamma2 vs alpha2/beta2/gamma2 median across exact Tier-1 Ki contexts for label candidate; raw context values retained candidate_only_not_primary_outcome
 triazolam        0.8       0.59                  0.737500                       -0.132238                          0.8                           0.59                                                                     0.7375                         1                                            alpha1/beta3/gamma2 vs alpha2/beta3/gamma2 median across exact Tier-1 Ki contexts for label candidate; raw context values retained candidate_only_not_primary_outcome
  zolpidem       27.0     156.00                  5.842697                        0.766613 27;26.7;52.1;26.7;27;27;26.7  160;156;255.1;156;103;160;156                       5.92593;5.8427;4.89635;5.8427;3.81481;5.92593;5.8427                         7                                            alpha1/beta3/gamma2 vs alpha2/beta3/gamma2 median across exact Tier-1 Ki contexts for label candidate; raw context values retained candidate_only_not_primary_outcome

Label medians aggregate exact Tier-1 Ki contexts only and retain raw context values.

## ML dataset readiness

                      status  n_drugs  n_features  missing_rate_all_features  zero_variance_feature_count  highly_sparse_frequency_feature_count  distance_stat_missing_cells  distance_stat_zero_cells  unresolved_interaction_rows_excluded  unresolved_feature_count  label_rows                                                         feature_leakage_risk                                                                                                                reason
NOT READY FOR MODEL TRAINING        4         200                    0.15625                           30                                      2                          100                         0                                   177                        17           4 low by column content; future ML must split by drug and keep labels out of X n=4, provisional alpha2 structure, unresolved mapping exclusions, and distance missingness require QC before training

NOT READY FOR MODEL TRAINING. The matrix is X-shaped, but n=4, provisional alpha2 structure, unresolved mapping exclusions, and distance missingness prevent model training or performance claims.

## Requested interpretation

1. Alpha1 fingerprint differences are visible in the primary alpha1 drug × feature matrix and pairwise table.
2. Alpha2 also shows drug-to-drug differences, conditional on the provisional 9CTJ-derived construct.
3. Highest-variance residue/interaction features are ranked in alpha1_variable_features.csv and alpha2_variable_features.csv; these are descriptive candidates.
4. Distance differences can be inspected where PLIP supplied a distance. Distance NA means no distance observation or no interaction, never zero distance.
5. Zolpidem high/low features are listed by drug_high and drug_low; uniqueness is not claimed.
6. Triazolam high/low features are listed in the same tables; uniqueness is not claimed.
7. Separate alpha1 and alpha2 blocks preserve receptor-context information that a difference-only representation would discard.
8. The combined matrix is structurally suitable as a future X table, but readiness is NOT READY FOR MODEL TRAINING.
9. Before adding drugs, resolve chain-aware mapping, review the provisional alpha2 model/box transfer, audit distance schemas, and keep beta2 pharmacology versus beta3 docking differences explicit.
10. The most natural next phase is to expand the drug panel with composition-matched pharmacology while preserving this feature schema; phenotype linkage should follow after feature/QC stability is demonstrated.

## QC boundaries

Pose counts: {('alprazolam', 'alpha1_beta3_gamma2'): 9, ('alprazolam', 'alpha2_beta3_gamma2_local_9CTJ'): 9, ('diazepam', 'alpha1_beta3_gamma2'): 9, ('diazepam', 'alpha2_beta3_gamma2_local_9CTJ'): 9, ('triazolam', 'alpha1_beta3_gamma2'): 9, ('triazolam', 'alpha2_beta3_gamma2_local_9CTJ'): 9, ('zolpidem', 'alpha1_beta3_gamma2'): 9, ('zolpidem', 'alpha2_beta3_gamma2_local_9CTJ'): 9}. All combinations currently have 9 poses. Frequency zero is used only for resolved features with pose data and no observed interaction. Unresolved mapping is excluded and listed in unresolved_feature_qc.csv. Distance missing is separate from distance=0. Alpha2 is provisional 9CTJ-derived; alprazolam pharmacology is beta2-background while docking is beta3-background. Existing strict and 4-drug reports/figures were not overwritten.

## Files

processed/alpha1_drug_feature_matrix.csv; alpha2_drug_feature_matrix.csv; combined_alpha1_alpha2_feature_matrix.csv; drug_pharmacology_labels.csv
tables/alpha1_interdrug_similarity.csv; alpha2_interdrug_similarity.csv; alpha1_variable_features.csv; alpha2_variable_features.csv; unresolved_feature_qc.csv; ml_dataset_readiness_qc.csv
figures/alpha1_interdrug_fingerprint_heatmap.png; alpha2_interdrug_fingerprint_heatmap.png; combined_multireceptor_fingerprint_heatmap.png; interaction_frequency_heatmap.png; interaction_distance_heatmap.png
