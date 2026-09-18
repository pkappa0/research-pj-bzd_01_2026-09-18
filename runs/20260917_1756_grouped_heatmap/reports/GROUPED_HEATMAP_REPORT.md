Run ID: 20260917_1756_grouped_heatmap
Parent Run: 20260917_1755_grouped_heatmap
Analysis Phase: residue-class grouped fingerprint visualization
Started: 2026-09-17T17:56:57+09:00
Completed: 2026-09-17T17:56:58+09:00

# GROUPED_HEATMAP_REPORT

The original residue_ligand_frequency_heatmap was not modified. This additional view uses the parent run's raw frequency matrix and changes only row ordering, labeling and the optional row mean-centered color scale.

## Residue-class visualization

Residues are classified for visualization only: Aromatic (PHE/TYR/TRP/HIS), Hydrophobic / aliphatic (ALA/VAL/LEU/ILE/MET/PRO), Hydrophilic / polar (SER/THR/ASN/GLN/CYS), Charged (LYS/ARG/ASP/GLU), and Other. Interaction type remains an independent feature dimension.

Rows are grouped into separate alpha1 and alpha2 blocks. Within each block, residue class, common position, and interaction type determine order; interaction types for the same residue are adjacent. The primary heatmap retains raw 0–1 frequencies.

## Row-normalized view

The normalized heatmap uses row mean-centering across the four drugs: each raw frequency minus that row's four-drug mean. It is a visualization of relative pattern only; the raw matrix is preserved unchanged and should be used for numeric reporting.

## Residue-specific pattern differences

            interpretation_group receptor common_position residue_name       residue_class        interaction_type  frequency_range highest_frequency_drug lowest_frequency_drug                                         drug_frequency_pattern                                                                   interpretation_note
      aromatic_pi_stack_features   alpha1  BZD_GAMMA2_058          TYR            Aromatic                pi_stack         0.333333              triazolam              zolpidem diazepam=0.444;alprazolam=0.444;triazolam=0.556;zolpidem=0.222 descriptive within-row pattern; no total-interaction ranking or statistical inference
      aromatic_pi_stack_features   alpha1  BZD_GAMMA2_077          PHE            Aromatic                pi_stack         0.333333             alprazolam             triazolam diazepam=0.222;alprazolam=0.444;triazolam=0.111;zolpidem=0.111 descriptive within-row pattern; no total-interaction ranking or statistical inference
      aromatic_pi_stack_features   alpha2    BZD_SITE_102          HIS            Aromatic                pi_stack         0.333333             alprazolam             triazolam diazepam=0.222;alprazolam=0.333;triazolam=0.000;zolpidem=0.000 descriptive within-row pattern; no total-interaction ranking or statistical inference
hydrophobic_interaction_features   alpha2    BZD_SITE_102          HIS            Aromatic hydrophobic_interaction         0.555556              triazolam              diazepam diazepam=0.000;alprazolam=0.111;triazolam=0.556;zolpidem=0.222 descriptive within-row pattern; no total-interaction ranking or statistical inference
hydrophobic_interaction_features   alpha1  BZD_GAMMA2_077          PHE            Aromatic hydrophobic_interaction         0.444444               zolpidem             triazolam diazepam=0.778;alprazolam=0.778;triazolam=0.444;zolpidem=0.889 descriptive within-row pattern; no total-interaction ranking or statistical inference
hydrophobic_interaction_features   alpha1    BZD_SITE_100          PHE            Aromatic hydrophobic_interaction         0.444444               diazepam             triazolam diazepam=0.444;alprazolam=0.333;triazolam=0.000;zolpidem=0.000 descriptive within-row pattern; no total-interaction ranking or statistical inference
   hydrophilic_or_hbond_features   alpha1  BZD_GAMMA2_060          ASN Hydrophilic / polar           hydrogen_bond         0.333333               zolpidem              diazepam diazepam=0.000;alprazolam=0.000;triazolam=0.111;zolpidem=0.333 descriptive within-row pattern; no total-interaction ranking or statistical inference
   hydrophilic_or_hbond_features   alpha1    BZD_SITE_205          SER Hydrophilic / polar           hydrogen_bond         0.222222               diazepam             triazolam diazepam=0.222;alprazolam=0.222;triazolam=0.000;zolpidem=0.000 descriptive within-row pattern; no total-interaction ranking or statistical inference
   hydrophilic_or_hbond_features   alpha2    BZD_SITE_205          SER Hydrophilic / polar           hydrogen_bond         0.111111               diazepam            alprazolam diazepam=0.111;alprazolam=0.000;triazolam=0.000;zolpidem=0.000 descriptive within-row pattern; no total-interaction ranking or statistical inference

These entries are selected by within-row frequency range. They are not ranked by total interaction count and do not imply statistical significance, subtype causality, or pharmacological potency.

## Key-site visibility

The row labels explicitly retain gamma2 TYR58, gamma2 PHE77, alpha-side HIS102, LYS156/157 positions where present, and SER205/206 positions where present. A row appears only when the corresponding feature was observed in the lossless parent matrix; no unobserved feature was added.

## Files

results/figures/residue_ligand_frequency_heatmap_grouped.png; results/figures/residue_ligand_frequency_heatmap_row_normalized.png; results/tables/residue_ligand_frequency_matrix_grouped.csv; results/tables/grouped_heatmap_top_differences.csv
