# Task: Create a wobbling-focused PLIF–clinical phenotype figure

Existing frozen structural fingerprint X and clinical fingerprint Y must be used as-is.
Do NOT recompute docking, PLIF, clinical ROR, feature definitions, or perform feature selection.

The purpose of this figure is to focus the existing exploratory analysis on the clinical phenotype "ふらつき / motor instability".

## Scientific framing

The observational unit is the drug.

We are NOT claiming:

    residue -> adverse event
    PLIF feature -> causation
    mediation through rotarod

We are asking:

    Across the same drugs,
    do variations in frozen PLIF features correspond to variations
    in clinically observed balance / motor-impairment phenotypes?

All analyses remain exploratory drug-level correspondence analyses.

---

# 1. Clinical phenotype grouping

Use the existing prespecified clinical terms and group them as follows.

## PRIMARY: Core balance / ataxia phenotype (C1)

- Ataxia
- Balance disorder
- Coordination abnormal
- Gait disturbance

This is the primary clinical phenotype for the present figure.

## SECONDARY: Motor weakness phenotype (C5)

- Muscular weakness
- Hypotonia

Treat C5 as a possible motor contributor to clinical wobbling,
not as identical to C1.

## CONTEXT ONLY

C3:
- Somnolence
- Sedation

C4:
- Dizziness
- Vertigo

C2:
- Fall

These are NOT part of the primary C1 endpoint.

Interpretation:

C3 = sedation-related contributor / confounder
C4 = subjective dizziness-related contributor / confounder
C2 = downstream clinical outcome

C6 cognitive terms should not be shown in the main wobbling-focused
correspondence heatmap.

---

# 2. Figure layout

Create one publication-quality figure:

    wobbling_focused_plif_clinical_correspondence.png
    wobbling_focused_plif_clinical_correspondence.svg

Use a multi-panel layout.

========================================
Panel A — Clinical phenotype definition
========================================

Create a simple conceptual diagram showing:

                    WOBBLING / MOTOR INSTABILITY

                 PRIMARY
        C1 Balance / Ataxia
        - Ataxia
        - Balance disorder
        - Coordination abnormal
        - Gait disturbance

                 SECONDARY
        C5 Motor weakness
        - Muscular weakness
        - Hypotonia

Then show smaller contextual boxes:

        C3 Sedation
        C4 Dizziness
        C2 Fall

Label them:

C3/C4:
"contributing / potentially confounding phenotypes"

C2:
"downstream outcome"

Do NOT draw causal arrows.

Use wording such as:

"Phenotype organization for exploratory correspondence analysis"

========================================
Panel B — PLIF vs wobbling phenotype
========================================

Using the existing archived/frozen drug-level Spearman results,
make a heatmap with:

Rows:
    all 48 frozen PLIF features F01-F48

Columns, in exactly this order:

    C1 Ataxia
    C1 Balance disorder
    C1 Coordination abnormal
    C1 Gait disturbance

    C5 Muscular weakness
    C5 Hypotonia

Separate C1 and C5 visually with a vertical gap or thick divider.

Cell color:
    archived Spearman rho

Color scale:
    -1 to +1
    diverging color map centered at 0

Cell text:
    n_effective

Gray:
    NA / undefined correlations

Do NOT rank features by rho.
Do NOT reorder features based on clinical results.
Keep frozen F01-F48 order.

Add a clear annotation:

"Each correlation is calculated across drugs.
Rows are structural PLIF features, not causal residues."

Use title:

"Drug-level PLIF correspondence with balance / motor phenotypes"

========================================
Panel C — Representative structural patterns
========================================

Because F02, F04, F16 and F28 appeared visually variable in the
frozen structural fingerprint, show their original drug-level PLIF
values as descriptive examples.

IMPORTANT:
These are NOT selected as statistically significant features.
Do not call them biomarkers or hits.

Make a small heatmap:

Rows:
    10 drugs in the same frozen order

Columns:
    F02
    F04
    F16
    F28

Use original frozen contact-frequency values.

Label:

"Representative visually variable PLIF features
(pre-specified for visualization only; no feature selection)"

Do NOT remove the other 44 features from Panel B.

========================================
Panel D — External in vivo anchor
========================================

Reuse the archived Nishino 2008 diazepam / triazolam rotarod data.

Plot:

diazepam vs triazolam

2 mg/kg:
15, 30, 60, 90 min

5 mg/kg:
15, 30, 60, 90 min

y-axis:
Rotarod failure fraction

Keep the data exactly as already archived.

Label the panel:

"External biological anchor: Nishino 2008"

Add:

"Context-matched diazepam/triazolam comparison only (n=2 drugs)"
"External anchor, not a mediator"
"No correlation or causal inference"

---

# 3. Optional small context panel

If space permits, add a narrow heatmap for:

C3 Somnolence
C3 Sedation
C4 Dizziness
C4 Vertigo
C2 Fall

Use all F01-F48 and the same archived Spearman values.

Clearly title:

"Context phenotypes — not primary wobbling endpoint"

Do NOT mix these columns into the primary C1/C5 interpretation.

---

# 4. Interpretation rules

Do not automatically identify individual PLIF features as mechanistic residues.

The figure should support statements only at this level:

"Across drugs, some PLIF features show concordant or discordant
variation with clinical balance/motor phenotype terms."

Do NOT write:

"Fxx causes ataxia"
"Fxx is the ataxia residue"
"PLIF predicts falls"
"Rotarod mediates the clinical phenotype"

If a PLIF feature shows similar rho direction across several C1 terms,
it may be described only as:

"a candidate feature showing internally consistent drug-level
correspondence across multiple C1 terms."

No p-value-based feature ranking.
No multiple-testing-driven hit calling in this PoC.

---

# 5. QC table

Also export:

    wobbling_plif_clinical_correspondence.csv

Columns:

feature_id
feature_label
clinical_axis
clinical_term
spearman_rho
n_effective

Add:

    c1_direction_consistency_summary.csv

For each F01-F48 calculate descriptive values only:

- number of estimable C1 terms
- number positive rho
- number negative rho
- median rho across C1 terms
- min rho
- max rho

Do NOT convert this into a significance score or feature ranking.

The purpose is to inspect whether individual structural features
show directionally consistent correspondence across the four
related C1 terms.

---

# 6. Reproducibility

Reuse existing frozen source files wherever possible.

Do not overwrite historical outputs.

Save the plotting / summary script in the repository, for example:

    scripts/plot_wobbling_focused_correspondence.py

Save outputs into the current analysis/results directory.

Document:
- exact input files
- column mappings
- drug order
- feature order
- missing-value handling
- output paths

Update the relevant README / analysis note with a short methods section.

Do not alter frozen X or Y.