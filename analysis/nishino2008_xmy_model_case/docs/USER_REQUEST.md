# MASTER TASK
# Nishino 2008 matched X–M–Y model case
# Structural → in vivo → Clinical
#
# Goal:
# Complete one biologically matched model case before expanding to ML.

============================================================
0. SCIENTIFIC PURPOSE
============================================================

The previous PoC established that:

Structural X:
- multi-receptor PLIF can represent drug-specific interaction patterns
- archived Vina docking scores provide an independent docking-score view

Clinical Y:
- adverse-event terms can be represented as a drug-level clinical fingerprint
- exploratory PLIF ↔ Clinical correspondence can be calculated using drug identity
  as the observational unit

However, in vivo M has remained sparse.

The present task should NOT attempt machine learning.

Instead, use Nishino et al. 2008 as a matched model case to determine whether
a coherent three-layer representation can be constructed:

    Structural X
         ↕
    in vivo M
         ↕
    Clinical Y

The primary research question is:

"Can structural, controlled in vivo, and clinical phenotype information
be aligned for the same drugs in a reproducible data representation?"

This is a model-case / representation-establishment study.

Do NOT try to prove causality.
Do NOT optimize Structural X against Clinical Y.
Do NOT perform supervised ML with n=4.

============================================================
1. PRIMARY DRUG SET
============================================================

Primary matched model-case drugs:

- diazepam
- triazolam
- brotizolam
- lormetazepam

These four are the primary set.

Nishino 2008 also includes rilmazafone.

DO NOT include rilmazafone in the primary structural model case because
rilmazafone is a prodrug and oral in vivo activity reflects active metabolites.

Create:

docs/rilmazafone_exclusion_note.md

Explain:

- Nishino 2008 contains rilmazafone
- it is retained in the raw in vivo source table
- it is excluded from the primary X–M–Y model case
- parent-compound docking should not be equated with the pharmacological
  activity produced after metabolic activation

Do not discard rilmazafone data.

It may later be used as a metabolism-aware sensitivity case, but NOT in
the primary analysis.

============================================================
2. SOURCE VERIFICATION — NISHINO 2008
============================================================

Primary paper:

Nishino et al. 2008
"Evaluation of Anxiolytic-like Effects of Some Short-Acting
Benzodiazepine Hypnotics in Mice"
Journal of Pharmacological Sciences
DOI: 10.1254/jphs.08107FP

Use the paper itself as the authoritative source.

Do NOT rely on values copied from previous chat summaries when the paper
can be checked directly.

Extract and document:

- species
- strain
- sex if reported
- administration route
- vehicle
- rotarod rotation speed
- trial duration
- pretraining / selection procedure
- observation times
- number of animals
- positive-event definition
- dose levels
- positive / tested counts
- reported ED50 values
- reported ED50 confidence intervals
- ED50 calculation method

Create:

data/nishino2008_rotarod_raw.csv

Required columns:

study_id
drug
dose_mg_kg
time_min
n_positive
n_total
failure_fraction
administration_route
species
strain
sex
rotarod_rpm
trial_duration_sec
source_table
source_reference

Retain ALL FIVE Nishino drugs in this raw source table.

Then create:

data/nishino2008_rotarod_ed50.csv

Columns:

drug
ed50_mg_kg
ed50_ci_low
ed50_ci_high
ed50_method
primary_model_case
exclusion_reason

Do not recompute the reported ED50 as the primary value.

If a derived ED50 is computed independently from Table 1,
store it separately and label it as a secondary validation only.

============================================================
3. DEFINE THE IN VIVO M LAYER
============================================================

For the primary four drugs, define:

M_raw = reported rotarod ED50 (mg/kg)

Lower ED50 = greater potency for rotarod impairment.

Also create a monotonic descriptive representation:

rotarod_potency = -log10(ED50_mg_kg)

This is ONLY for intuitive visualization where:

higher value = greater rotarod impairment potency

Keep the original ED50 as the primary source value.

Do not call rotarod ED50:

- clinical risk
- ataxia incidence
- binding affinity

Use:

"rotarod motor-impairment potency"
or
"rotarod phenotype"

Also retain the complete dose × time response matrix separately.

Create:

data/nishino2008_primary4_in_vivo_M.csv

============================================================
4. STRUCTURAL X — CHECK EXISTING STANDARDIZED ENVIRONMENT
============================================================

Existing standardized structural pipeline contains:

- alpha1 receptor block
- alpha2 receptor block
- identical docking box
- identical docking parameters
- 5 independent seeds
- up to 9 poses per seed
- PLIF generation
- archived Vina score extraction

Existing standardized 10-drug run already contains:

- diazepam
- triazolam

It does NOT currently contain:

- brotizolam
- lormetazepam

Before adding the new drugs, verify that the exact standardized environment
can be reproduced.

Check and record:

- receptor files / hashes
- docking box
- exhaustiveness
- num_modes
- seed handling
- Vina version
- ligand-preparation procedure
- receptor-preparation procedure
- PLIP version
- PLIF extraction logic
- residue mapping
- software environment / container if available

Create:

qc/nishino_structural_environment_check.md

DECISION RULE:

If the standardized environment is fully reproducible:
    reuse existing frozen diazepam and triazolam outputs
    and add only brotizolam and lormetazepam.

If there is ANY material uncertainty that the new two drugs can be generated
under exactly the same environment:
    create a new isolated "Nishino model-case" run and redock ALL FOUR drugs.

Do not mix incompatible docking batches.

Document the decision explicitly.

============================================================
5. PREPARE BROTIZOLAM AND LORMETAZEPAM
============================================================

For brotizolam and lormetazepam:

- obtain canonical chemical identity from an authoritative chemistry source
- record identifiers and canonical SMILES
- generate ligand preparation using the exact standardized procedure
- verify protonation / charge handling used by the existing pipeline
- do not manually optimize structures based on clinical outcomes

Create:

data/nishino_primary4_ligand_manifest.csv

Columns:

drug
canonical_smiles
source_database
source_identifier
structure_preparation_method
charge_method
notes

============================================================
6. STANDARDIZED DOCKING
============================================================

Generate / reuse structural data for:

4 drugs
× 2 receptor blocks
× 5 independent seeds
× up to 9 poses

Primary receptors:

- alpha1
- alpha2

Do NOT add alpha3 or alpha5 in this model-case task.

For every output retain:

drug
receptor
seed
pose_rank
vina_score
source_pdbqt

Run the same QC previously used in the standardized 10-drug analysis.

============================================================
7. PLIF GENERATION
============================================================

Generate the same frozen PLIF representation.

Do NOT redesign the fingerprint after seeing Nishino ED50 or Clinical Y.

Do NOT select residues based on rotarod or FAERS results.

Use the existing residue mapping / feature definition.

The structural representation remains:

X_PLIF =
[
    alpha1 PLIF
    |
    alpha2 PLIF
]

Create the four-drug matrix:

data/nishino_primary4_plif.csv

Preserve the same F01-F48 feature identity wherever technically valid.

If new docking outputs produce residue-feature incompatibilities,
do not silently modify the fingerprint.

Flag them.

============================================================
8. VINA DOCKING-SCORE LAYER
============================================================

Parse / reuse all archived:

REMARK VINA RESULT

records.

For each drug × receptor × seed:

identify the best-scoring pose.

For each drug × receptor calculate:

- n_seeds
- median best-seed Vina score
- mean
- SD
- IQR
- min
- max

Primary structural score:

median_best_seed_score

Create:

data/nishino_primary4_vina_summary.csv

Also create:

vina_favorability = -median_best_seed_score

Retain the raw score.

Terminology:

"Vina docking score"

NOT:

"experimental affinity"
"true affinity"
"binding free energy"

============================================================
9. STRUCTURAL X DEFINITION
============================================================

The primary structural representation should remain multi-view:

X_structural =
{
    PLIF_alpha1,
    PLIF_alpha2,
    Vina_alpha1,
    Vina_alpha2
}

Do NOT collapse these into a single composite score.

Do NOT multiply PLIF by Vina score.

Do NOT weight structural features using Clinical Y.

Create:

data/nishino_primary4_structural_multiview.csv

One row per drug.

============================================================
10. CLINICAL Y — COVERAGE AUDIT FIRST
============================================================

Use the EXACT existing frozen Clinical-fingerprint methodology.

Do not redefine adverse-event terms because of the Nishino drugs.

Primary raw C1 candidate terms remain:

- Ataxia
- Balance disorder
- Coordination abnormal
- Gait disturbance

The C1 grouping is PROVISIONAL.

The MedDRA PT-level values are the actual raw clinical variables.

For diazepam and triazolam:
reuse existing frozen Y if applicable.

For brotizolam and lormetazepam:
apply the exact same clinical-data extraction pipeline.

Before interpreting Y, create a coverage audit.

Create:

qc/nishino_primary4_clinical_coverage.csv

For each:

drug × MedDRA PT

record:

- exposed / target report count used by the existing ROR pipeline
- comparator count(s)
- total reports
- estimable / not estimable
- missingness reason
- logROR if estimable
- exact drug normalization / synonym rule used

IMPORTANT:

Do NOT impute missing clinical signals.

Do NOT substitute zero for missing data.

Do NOT silently switch databases for individual drugs.

If brotizolam or lormetazepam has insufficient FAERS/openFDA coverage,
state this explicitly.

The model case is allowed to produce:

complete X–M
but incomplete X–M–Y

if that is what the data support.

Do not manufacture a full three-layer result.

============================================================
11. CLINICAL Y TABLE
============================================================

If coverage is adequate, create:

data/nishino_primary4_clinical_Y.csv

Primary variables:

drug
logROR_ataxia
logROR_balance_disorder
logROR_coordination_abnormal
logROR_gait_disturbance

Optionally preserve related terms as SECONDARY CONTEXT:

- muscular weakness
- hypotonia
- dizziness
- vertigo
- somnolence
- sedation
- fall

Do not combine them into a single score in the primary analysis.

============================================================
12. BUILD THE MATCHED X–M–Y TABLE
============================================================

Create the central model-case table:

data/nishino_primary4_XMY_master.csv

One row per drug.

Required conceptual blocks:

IDENTITY
- drug

STRUCTURAL X
- alpha1 PLIF reference
- alpha2 PLIF reference
- alpha1 median Vina score
- alpha2 median Vina score
- alpha1 Vina favorability
- alpha2 Vina favorability

IN VIVO M
- rotarod ED50
- ED50 CI
- rotarod potency

CLINICAL Y
- Ataxia logROR
- Balance disorder logROR
- Coordination abnormal logROR
- Gait disturbance logROR

QC
- structural completeness
- M completeness
- Y completeness

============================================================
13. PRIMARY FIGURE — MATCHED MODEL CASE
============================================================

Create:

figures/nishino_primary4_XMY_model_case.png
figures/nishino_primary4_XMY_model_case.svg

Rows:

diazepam
triazolam
brotizolam
lormetazepam

Keep this exact study-defined order or document another fixed order.

Create horizontally aligned blocks:

A. Structural PLIF
   - alpha1
   - alpha2
   - F01-F48 frozen order
   - contact-frequency heatmap

B. Vina docking score
   - alpha1
   - alpha2
   - raw score annotated
   - color may use within-receptor favorability/rank
   - do not imply raw alpha1 and alpha2 scores are directly comparable

C. in vivo M
   - reported rotarod ED50
   - derived rotarod potency
   - show ED50 CI if visually possible

D. Clinical Y
   - Ataxia
   - Balance disorder
   - Coordination abnormal
   - Gait disturbance
   - logROR heatmap
   - gray = unavailable

Title:

"Nishino 2008 matched model case:
multi-view structural, in vivo, and clinical phenotype representation"

Subtitle:

"Four direct-acting benzodiazepines aligned by drug identity"

Footnote:

"Representation / hypothesis-generation only.
No causal, mediation, or predictive inference."

============================================================
14. FIGURE — ROTAROD PHENOTYPE
============================================================

Create a dedicated Nishino rotarod figure.

Use all raw Nishino 2008 Table 1 values for the primary four drugs.

Show dose × time behavior.

Preferred presentation:

one small panel per drug

x = time (15, 30, 60, 90 min)
y = failure fraction
line = dose

OR a clear dose × time heatmap.

Also show reported ED50 separately.

Create:

figures/nishino_primary4_rotarod_profiles.png
figures/nishino_primary4_rotarod_profiles.svg

Do not normalize doses across drugs in the raw panel.

============================================================
15. DESCRIPTIVE X ↔ M ANALYSIS
============================================================

Because n=4:

NO supervised ML.
NO significance-based feature discovery.
NO biomarker calling.

Perform descriptive exploratory correspondence only.

For Vina:

calculate Spearman across the four drugs:

alpha1 Vina favorability ↔ rotarod potency
alpha2 Vina favorability ↔ rotarod potency

For each PLIF feature:

calculate descriptive Spearman:

PLIF feature ↔ rotarod potency

BUT:

- report n=4 prominently
- do not rank features as discoveries
- do not use p-value hit calling
- do not alter the fingerprint based on these values

Create:

data/nishino_primary4_X_M_correspondence.csv

and a compact descriptive figure:

figures/nishino_primary4_X_M_correspondence.png

Use wording:

"descriptive correspondence in a four-drug matched model case"

============================================================
16. DESCRIPTIVE M ↔ Y ANALYSIS
============================================================

If Clinical Y is sufficiently complete:

calculate:

rotarod potency ↔ each C1 clinical PT

Across the primary four drugs.

Again:

n <= 4
descriptive only.

Create:

data/nishino_primary4_M_Y_correspondence.csv

Do NOT call this validation of causality.

The question is only:

"Does the ordering of controlled motor-impairment potency show any
directional agreement with clinical reporting signals?"

If coverage is inadequate, produce the table with NA and explain why.

============================================================
17. DESCRIPTIVE X ↔ Y ANALYSIS
============================================================

Reuse the same logic already established in the 10-drug PoC.

Compare:

PLIF ↔ C1 clinical terms
Vina favorability ↔ C1 clinical terms

for the four-drug model case.

This is SECONDARY.

Do not replace the existing 10-drug X ↔ Y result.

The 10-drug analysis remains the stronger clinical correspondence PoC.

The four-drug analysis exists only to connect X–M–Y within one matched
in vivo study.

============================================================
18. DO NOT CREATE A MEDIATION MODEL
============================================================

Do NOT run:

X -> M -> Y mediation

Do NOT claim:

- M mediates X and Y
- a structural feature causes rotarod impairment
- rotarod predicts clinical fall/ataxia
- one receptor subtype is dominant
- a residue is protective
- a residue causes wobbling

With n=4, none of these are justified.

Use:

"matched three-layer representation"
"directional correspondence"
"model-case"
"working hypothesis"

============================================================
19. FINAL SYNTHESIS FIGURE
============================================================

Create one clean conceptual + empirical summary figure:

figures/nishino_three_layer_summary.png
figures/nishino_three_layer_summary.svg

Concept:

               Structural X
     [PLIF α1 | PLIF α2 | Vina α1 | Vina α2]
                       ↕
                 in vivo M
          [Nishino rotarod ED50]
                       ↕
                 Clinical Y
     [Balance / Ataxia-related MedDRA PTs]

Alongside the concept, show whether each layer is:

- complete
- partially complete
- unavailable

for each of the four drugs.

Do NOT use causal arrows.

============================================================
20. MODEL-CASE SUCCESS CRITERIA
============================================================

Evaluate the model case using PREDEFINED process-level criteria.

Do NOT define success as obtaining a positive correlation.

Evaluate:

A. Structural feasibility
Can a standardized multi-view X be created for all four drugs?

B. in vivo feasibility
Can a common controlled M parameter be assigned to all four drugs
from one study?

C. Clinical feasibility
Can comparable Clinical Y terms be estimated for the same drugs?

D. Alignment feasibility
Can all available layers be joined by unambiguous drug identity?

E. Biological dynamic range
Do the layers contain sufficient drug-to-drug variability to justify
a larger study?

F. No circularity
Was X constructed independently of M and Y?

The model case is considered informative even if correlations are weak.

============================================================
21. GO / NO-GO DECISION FOR NEXT STAGE
============================================================

Create:

docs/nishino_model_case_go_nogo.md

Do not give a simplistic PASS / FAIL.

Instead evaluate three possible outcomes:

CASE 1:
X, M, and Y are all sufficiently populated.

Conclusion:
Proceed to systematic expansion of matched in vivo datasets.

CASE 2:
X and M are strong, but Y is sparse for brotizolam / lormetazepam.

Conclusion:
The bottleneck is clinical-data coverage.
Do NOT interpret this as a biological failure.
Next milestone is selection / validation of an appropriate clinical
data source before expansion.

CASE 3:
Structural or in vivo representation itself cannot be harmonized.

Conclusion:
Resolve representation / assay-context problems before scaling.

============================================================
22. FUTURE IN VIVO EXPANSION PLAN
============================================================

Create:

docs/future_in_vivo_expansion_plan.md

Based on what was learned from the Nishino model case, define the minimum
metadata required for future literature-derived in vivo datasets:

- drug identity
- species
- strain
- sex
- administration route
- dose
- time after administration
- assay protocol
- rotarod speed / acceleration
- training procedure
- endpoint definition
- n
- effect value
- uncertainty / CI where available

Separate:

Tier 1:
same study, multiple drugs, same experimental context

Tier 2:
closely matched protocols across studies

Tier 3:
heterogeneous studies unsuitable for direct numeric comparison

Do NOT pool Tier 3 values as if equivalent.

============================================================
23. FUTURE ML ROADMAP
============================================================

Create:

docs/future_ml_roadmap.md

The current task must NOT implement ML.

Instead define the future sequence:

Stage 0 — current
Establish X / M / Y representation and data schema.

Stage 1
Systematically collect matched in vivo datasets.

Stage 2
Expand drug n while keeping raw Clinical PT-level outcomes.

Stage 3
Expand receptor representation only after sufficient drug n exists
(e.g. alpha3 / alpha5).

Stage 4
Compare predictive baselines:

Model A:
Vina only -> Clinical Y

Model B:
PLIF only -> Clinical Y

Model C:
PLIF + Vina -> Clinical Y

Model D:
PLIF + Vina + available in vivo M -> Clinical Y

Stage 5
Investigate methods for sparse M:

- multi-task learning
- multi-view learning
- missing-label / partially observed intermediate phenotype approaches
- latent representations

Do not select an algorithm prematurely.

Stage 6
External drug validation.

Stage 7
If justified, use model uncertainty / information gain to prioritize
which additional drug should receive new in vivo characterization.

This is where active-learning / experimental-design logic may enter.

============================================================
24. FINAL REPORT
============================================================

Create:

FINAL_REPORT_NISHINO_XMY_MODEL_CASE.md

Required structure:

# 1. Research question

# 2. Why Nishino 2008 was selected

# 3. Drug-set definition
- primary four drugs
- rilmazafone exclusion rationale

# 4. Structural X
- PLIF
- Vina
- QC

# 5. in vivo M
- rotarod protocol
- raw dose/time data
- ED50
- limitations

# 6. Clinical Y
- source
- MedDRA PTs
- coverage
- missingness

# 7. Matched X–M–Y representation

# 8. Descriptive correspondence
- X ↔ M
- M ↔ Y
- X ↔ Y
- no causal interpretation

# 9. What this model case establishes

# 10. What it does NOT establish

# 11. Main bottleneck identified

# 12. Go / no-go reasoning

# 13. Next experimental / data milestone

# 14. Future ML role

# 15. Reproducibility
- exact inputs
- scripts
- versions
- outputs
- commit hash

The final conclusion must NOT depend on obtaining positive correlations.

The primary output should answer:

"Is a matched Structural–in vivo–Clinical data model technically and
scientifically constructible, and what is the limiting layer?"

============================================================
25. REQUIRED FINAL OUTPUT INDEX
============================================================

At the end, create:

OUTPUT_INDEX.md

List every generated:

- source table
- QC table
- structural table
- in vivo table
- clinical table
- matched XMY table
- figure
- script
- documentation file
- final report

For each output provide one-line purpose.

============================================================
26. REPRODUCIBILITY / GIT
============================================================

Do not overwrite historical 9/20 PoC results.

Create a dedicated analysis directory such as:

analysis/nishino2008_xmy_model_case/

Keep:

data/
figures/
qc/
scripts/
docs/

Commit all reproducible scripts and text outputs.

Large intermediate docking files may follow the repository's existing
data-management convention.

Record the final Git commit hash in FINAL_REPORT_NISHINO_XMY_MODEL_CASE.md.

============================================================
27. ABSOLUTE GUARDRAILS
============================================================

Do NOT:

- tune X against Y
- call FAERS logROR incidence
- call ROR a risk ratio
- call Vina score experimental affinity
- treat rilmazafone parent docking as equivalent to oral rilmazafone activity
- impute missing Clinical Y
- claim mediation
- claim causality
- perform Random Forest / neural network / supervised ML with n=4
- select "best residues" from the four-drug analysis
- interpret a positive correlation as model validation
- suppress negative results

The purpose is a fair matched-model-case evaluation.

Positive and negative results are equally valid.