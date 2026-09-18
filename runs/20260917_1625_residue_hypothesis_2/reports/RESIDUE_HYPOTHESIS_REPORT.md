Run ID: 20260917_1625_residue_hypothesis_2
Parent Run: 20260917_1604_multireceptor_feature
Analysis Phase: BZD interface residue hypothesis validation
Started: 2026-09-17T16:25:39+09:00
Completed: 2026-09-17T16:25:39+09:00

# RESIDUE_HYPOTHESIS_REPORT

This run tests residue-level, within-receptor interaction hypotheses. It does not perform global alpha1-vs-alpha2 similarity analysis, ML, imputation, activity integration, or docking.

## Recovery and mapping

Raw PLIP rows: **339**; previously unresolved rows: **177**; recovered through gamma2 C→E mapping: **177**; unresolved after correction: **0**.
Gamma2 alignment score: **343.0**; optimal alignments: **1**. The selected alignment is unique in the Biopython PairwiseAligner result.
The prior alpha1 chain D→alpha2 chain D mapping is retained. The new BZD-interface mapping is gamma2 chain C in 6HUP to gamma2 chain E in the provisional 9CTJ local construct. Mapping is sequence/structure based; no PDB-number-only identity was assumed.

## Feature policy

`hydrophobic_interaction` and `pi_stack` are separate features. Pi stacking retains centdist, angle, offset, and PLIP type P/T frequencies. Hydrophobic contacts retain PLIP distance statistics. Hydrogen bonds, halogen bonds, and other recovered interaction types remain in the same residue-level tables.

## TYR58/PHE77 aromatic interaction summary

receptor    drug_id   site_label        interaction_type  frequency  distance_mean  centdist_mean  angle_mean  offset_mean  type_P_frequency  type_T_frequency
  alpha1   diazepam gamma2_TYR58 hydrophobic_interaction   0.777778       3.516000            NaN         NaN          NaN          0.000000          0.000000
  alpha1   diazepam gamma2_PHE77 hydrophobic_interaction   0.777778       3.490000            NaN         NaN          NaN          0.000000          0.000000
  alpha1   diazepam gamma2_TYR58                pi_stack   0.444444            NaN         3.9675     18.6875       1.5375          0.444444          0.000000
  alpha1   diazepam gamma2_PHE77                pi_stack   0.222222            NaN         4.5250     17.1850       1.4000          0.222222          0.000000
  alpha1 alprazolam gamma2_TYR58 hydrophobic_interaction   1.000000       3.525333            NaN         NaN          NaN          0.000000          0.000000
  alpha1 alprazolam gamma2_PHE77 hydrophobic_interaction   0.777778       3.502500            NaN         NaN          NaN          0.000000          0.000000
  alpha1 alprazolam gamma2_TYR58                pi_stack   0.444444            NaN         4.1775     16.5225       1.5325          0.444444          0.000000
  alpha1 alprazolam gamma2_PHE77                pi_stack   0.444444            NaN         4.5700     34.4550       1.3650          0.333333          0.111111
  alpha1  triazolam gamma2_TYR58 hydrophobic_interaction   0.888889       3.716667            NaN         NaN          NaN          0.000000          0.000000
  alpha1  triazolam gamma2_PHE77 hydrophobic_interaction   0.444444       3.693333            NaN         NaN          NaN          0.000000          0.000000
  alpha1  triazolam gamma2_TYR58                pi_stack   0.555556            NaN         4.1420     15.2800       1.6140          0.555556          0.000000
  alpha1  triazolam gamma2_PHE77                pi_stack   0.111111            NaN         4.1900     29.7800       1.3300          0.111111          0.000000
  alpha1   zolpidem gamma2_TYR58 hydrophobic_interaction   0.666667       3.734286            NaN         NaN          NaN          0.000000          0.000000
  alpha1   zolpidem gamma2_PHE77 hydrophobic_interaction   0.888889       3.562000            NaN         NaN          NaN          0.000000          0.000000
  alpha1   zolpidem gamma2_TYR58                pi_stack   0.222222            NaN         3.9600     15.0300       1.2050          0.222222          0.000000
  alpha1   zolpidem gamma2_PHE77                pi_stack   0.111111            NaN         4.7000     70.6000       0.4700          0.000000          0.111111
  alpha2   diazepam gamma2_TYR58 hydrophobic_interaction   0.333333       3.690000            NaN         NaN          NaN          0.000000          0.000000
  alpha2   diazepam gamma2_PHE77 hydrophobic_interaction   0.666667       3.577500            NaN         NaN          NaN          0.000000          0.000000
  alpha2   diazepam gamma2_TYR58                pi_stack   0.000000            NaN            NaN         NaN          NaN          0.000000          0.000000
  alpha2 alprazolam gamma2_TYR58 hydrophobic_interaction   0.444444       3.742500            NaN         NaN          NaN          0.000000          0.000000
  alpha2 alprazolam gamma2_PHE77 hydrophobic_interaction   0.666667       3.547143            NaN         NaN          NaN          0.000000          0.000000
  alpha2 alprazolam gamma2_TYR58                pi_stack   0.000000            NaN            NaN         NaN          NaN          0.000000          0.000000
  alpha2  triazolam gamma2_TYR58 hydrophobic_interaction   0.222222       3.860000            NaN         NaN          NaN          0.000000          0.000000
  alpha2  triazolam gamma2_PHE77 hydrophobic_interaction   0.444444       3.747500            NaN         NaN          NaN          0.000000          0.000000
  alpha2  triazolam gamma2_TYR58                pi_stack   0.111111            NaN         5.3900     88.0800       1.1100          0.000000          0.111111
  alpha2   zolpidem gamma2_TYR58 hydrophobic_interaction   0.444444       3.712500            NaN         NaN          NaN          0.000000          0.000000
  alpha2   zolpidem gamma2_PHE77 hydrophobic_interaction   0.666667       3.422500            NaN         NaN          NaN          0.000000          0.000000
  alpha2   zolpidem gamma2_TYR58                pi_stack   0.111111            NaN         5.3500     77.5300       1.5300          0.000000          0.111111

The priority-site rows above are descriptive within each receptor. A missing row for a drug is represented as `no_interaction_observed` when pose data exist; it is not treated as missing data or as a zero geometry.

## Drug-level residue comparison

receptor                    receptor_id common_position interface_subunit residue_name interaction_type  n_drugs_with_interaction  frequency_range frequency_high_drug frequency_low_drug                                                           frequency_values  distance_mean_range  centdist_mean_range  angle_mean_range  offset_mean_range                                                                                                                   profile_note
  alpha1            alpha1_beta3_gamma2  BZD_GAMMA2_058            gamma2          TYR         pi_stack                         4         0.333333           triazolam           zolpidem diazepam=0.444444;alprazolam=0.444444;triazolam=0.555556;zolpidem=0.222222                  NaN               0.2175            3.6575              0.409 frequency=interacting poses / available poses; geometry statistics use raw PLIP rows; hydrophobic and pi_stack remain separate
  alpha1            alpha1_beta3_gamma2  BZD_GAMMA2_077            gamma2          PHE         pi_stack                         4         0.333333          alprazolam          triazolam diazepam=0.222222;alprazolam=0.444444;triazolam=0.111111;zolpidem=0.111111                  NaN               0.5100           53.4150              0.930 frequency=interacting poses / available poses; geometry statistics use raw PLIP rows; hydrophobic and pi_stack remain separate
  alpha2 alpha2_beta3_gamma2_local_9CTJ  BZD_GAMMA2_058            gamma2          TYR         pi_stack                         2         0.000000           triazolam          triazolam               diazepam=0;alprazolam=0;triazolam=0.111111;zolpidem=0.111111                  NaN               0.0400           10.5500              0.420 frequency=interacting poses / available poses; geometry statistics use raw PLIP rows; hydrophobic and pi_stack remain separate

Frequency ranges compare the four drugs on the same receptor and feature. Geometry ranges are descriptive and are not used to infer a statistical effect.

## QC boundaries

Feature rows: alpha1=72, alpha2=76; corrected mapping rows: 302. Recovered interaction rows are retained in tables/recovered_interactions_qc.csv with source_row_id and raw_attributes.
The alpha2 receptor remains a provisional local construct derived from 9CTJ. Therefore alpha1-vs-alpha2 comparisons are secondary context only; the primary question here is drug-to-drug variation within alpha1 and within alpha2.
No interaction was collapsed across interaction types, residues, drugs, poses, or geometry fields. No values were imputed.

## Files

tables/corrected_residue_mapping.csv; tables/recovered_interactions_qc.csv; tables/residue_interaction_features_alpha1.csv; tables/residue_interaction_features_alpha2.csv; tables/residue_level_drug_comparison.csv; tables/aromatic_interaction_summary.csv
