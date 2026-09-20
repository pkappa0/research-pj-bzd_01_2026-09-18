# Task: Recover and analyze the archived Vina docking-score layer

## Purpose

The current frozen structural representation contains:

X_PLIF = [alpha1 PLIF | alpha2 PLIF]

We now want to recover the Vina docking scores already contained in the
archived standardized docking PDBQT outputs and establish an independent
"Affinity / docking-score layer".

This is NOT a new docking run.

Do NOT rerun docking.
Do NOT modify the frozen PLIF fingerprint.
Do NOT modify the frozen Clinical fingerprint.
Do NOT use Clinical Y to optimize, select, rescale, or tune structural features.

The goal is:

1. recover archived Vina scores reproducibly
2. QC score robustness across docking seeds
3. summarize alpha1 / alpha2 docking scores per drug
4. explore docking-score correspondence with the existing clinical phenotype
5. prepare the data for later PLIF + docking-score integration

---

# 1. Locate archived standardized docking outputs

Find the PDBQT docking outputs corresponding to the frozen
10-drug standardized redocking analysis.

Drugs:

- diazepam
- alprazolam
- triazolam
- zolpidem
- lorazepam
- clonazepam
- midazolam
- temazepam
- zopiclone
- zaleplon

Receptor blocks:

- alpha1
- alpha2

Expected design:

- 5 independent docking seeds per drug/receptor
- up to 9 poses per seed

Use ONLY outputs belonging to the standardized frozen 10-drug run.

Do not use older historical 4-drug PoC outputs.

If directory identity is ambiguous, inspect repository documentation,
run manifests, or filenames and document how the standardized run was identified.

---

# 2. Parse all Vina score records

From each docking PDBQT, parse all lines corresponding to:

REMARK VINA RESULT:

For every pose, retain:

- drug
- receptor_block
- seed
- pose_rank
- vina_score_kcal_mol
- rmsd_lb
- rmsd_ub
- source_pdbqt
- source_run / directory if useful

Save:

data/vina_scores_all_poses.csv

Do not discard non-best poses.

Flag:

- missing files
- malformed RESULT records
- missing seeds
- unexpected pose counts
- duplicated outputs

---

# 3. Seed-level best-score table

For each:

drug × receptor × seed

identify the best-scoring pose
(most negative Vina score).

Save:

data/vina_scores_seed_best.csv

Columns should include at minimum:

drug
receptor_block
seed
best_pose_rank
best_vina_score_kcal_mol
source_pdbqt

---

# 4. Drug × receptor summary

For each:

drug × receptor

summarize the five seed-level best scores.

Calculate:

- n_seeds_available
- median_best_seed_score
- mean_best_seed_score
- sd_best_seed_score
- IQR_best_seed_score
- min_best_seed_score
- max_best_seed_score

Primary descriptive value:

median_best_seed_score

Do NOT use the single globally best pose as the primary drug value.

Save:

data/vina_scores_drug_receptor_summary.csv

---

# 5. Important terminology

Use:

"Vina docking score"

Do NOT call it:

- measured binding affinity
- experimental affinity
- true binding free energy

More negative Vina score means a more favorable docking score.

For plots where intuitive direction would help, it is acceptable to add
a clearly derived column:

vina_favorability = -1 × median_best_seed_score

but retain the original raw Vina score everywhere.

Never replace or overwrite the raw value.

---

# 6. Figure A — Drug × receptor docking-score overview

Create:

figures/vina_score_alpha1_alpha2_heatmap.png
figures/vina_score_alpha1_alpha2_heatmap.svg

Rows:
10 drugs in the existing frozen drug order.

Columns:
alpha1
alpha2

Cell value:
median_best_seed_score

Annotate each cell with the numerical score.

Also indicate n_seeds if any value has fewer than 5 available seeds.

Title:

"Archived Vina docking scores across α1 and α2 receptor blocks"

Subtitle:

"Median of independent seed-level best scores; no redocking"

Do not merge alpha1 and alpha2.

---

# 7. Figure B — Seed robustness

Create a QC plot showing all five seed-level best scores for each
drug × receptor.

Preferred layout:

- x = drug
- y = Vina docking score
- separate alpha1 / alpha2 panels or otherwise clearly distinguish blocks
- show individual seed values
- optionally show median

Purpose:

Inspect whether a drug/receptor score is stable across independent seeds.

Save:

figures/vina_score_seed_robustness.png
figures/vina_score_seed_robustness.svg

---

# 8. Alpha1 vs alpha2 descriptive comparison

Create a scatter plot:

x = alpha1 median best-seed Vina score
y = alpha2 median best-seed Vina score

One point per drug.
Label each drug.

Add Spearman rho descriptively.

Do not interpret this as receptor equivalence or experimental affinity.

Save:

figures/vina_score_alpha1_vs_alpha2.png
figures/vina_score_alpha1_vs_alpha2.svg

Purpose:

Determine whether drugs that score favorably within alpha1
also tend to score favorably within alpha2.

---

# 9. Clinical correspondence — affinity layer only

Use the already frozen Clinical fingerprint Y.

Do NOT regenerate FAERS data.

Primary phenotype remains the current provisional C1 terms:

- Ataxia
- Balance disorder
- Coordination abnormal
- Gait disturbance

For each receptor independently:

calculate across drugs:

Spearman(
    median_best_seed_vina_score,
    clinical logROR
)

for each C1 term.

Because raw Vina score becomes MORE favorable as it becomes more negative,
also provide a clearly interpretable derived version using:

vina_favorability = -vina_score

The two contain identical rank information with reversed rho sign.

Use the favorability form for the main visualization,
while preserving raw-score correlations in the output table.

Save:

data/vina_clinical_c1_correspondence.csv

Include:

receptor_block
clinical_term
n_effective
rho_raw_vina_score
rho_vina_favorability

No p-value-based hit calling.
No causal interpretation.

---

# 10. Figure C — Affinity vs C1 clinical phenotype

Create a compact figure showing:

alpha1 Vina favorability ↔ C1 terms
alpha2 Vina favorability ↔ C1 terms

Suggested design:

rows:
alpha1
alpha2

columns:
Ataxia
Balance disorder
Coordination abnormal
Gait disturbance

cell color:
Spearman rho using vina_favorability

cell text:
n_effective

Title:

"Vina docking-score correspondence with C1 Balance / Ataxia phenotype"

Footnote:

"Observational unit = drug.
Vina score is a docking-score proxy, not experimental affinity.
Association/correspondence only."

Save PNG and SVG.

---

# 11. Do NOT integrate PLIF and Vina yet

For this task, keep:

PLIF layer
and
Vina docking-score layer

analytically separate.

Do NOT:

- weight PLIF by Vina score
- create a composite structural score
- optimize a PLIF+Vina metric against Clinical Y
- select residues based on affinity
- rerank PLIF features using clinical outcomes

The purpose of this task is to establish the independent affinity layer first.

---

# 12. Prepare for the next integration step

Create one machine-readable table:

data/structural_multiview_index.csv

One row per drug.

Include references / columns for:

- alpha1 PLIF vector identifier
- alpha2 PLIF vector identifier
- alpha1 median Vina score
- alpha2 median Vina score
- alpha1 Vina favorability
- alpha2 Vina favorability

Do not combine these numerically yet.

This will later support:

X_structural =
[
 alpha1 PLIF,
 alpha2 PLIF,
 alpha1 Vina,
 alpha2 Vina
]

---

# 13. Documentation

Update the relevant README / analysis note with:

- exact PDBQT directories used
- why they correspond to the standardized frozen run
- parser logic
- number of files / poses recovered
- any missing data
- seed summary definition
- distinction between Vina docking score and experimental affinity
- output files and figures

Also save the parsing / analysis script, e.g.:

scripts/recover_vina_score_layer.py

Do not overwrite historical outputs.
Commit the new outputs and documentation to the repository.