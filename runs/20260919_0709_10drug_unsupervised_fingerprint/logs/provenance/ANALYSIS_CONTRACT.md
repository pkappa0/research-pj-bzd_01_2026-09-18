# Analysis contract

## Primary comparison

**same receptor → different drugs**

Fix the receptor block and structural preparation, then compare ligands.
For the two-drug PoC, every primary drug difference is **alprazolam − diazepam**.
Report α1 and α2 blocks separately. Neither block is a replicate of the other.

## Secondary comparison

**same drug → different receptors**

Cross-receptor differences are supplementary and require explicit construct,
species, composition and preparation caveats. They must not replace the primary
within-receptor drug comparison or become its default plot axis.

## Multi-receptor representation

**drug = [α1 fingerprint | α2 fingerprint | ...]**

Use a fixed, versioned feature order and receptor-qualified feature names.
Concatenation preserves block identity; it does not average receptor blocks,
normalize away their identities, or establish a predictive model.
Retain interaction type and geometry with units. Preserve missingness masks.
Zero contact frequency requires an audited denominator and resolved mapping;
missing geometry is never zero. Geometry is conditional on observed contacts.

## Pharmacology endpoints

**binding / potency / efficacy are separate endpoints.**

Prefer different drugs in the same composition, species, assay, document,
endpoint, unit and relation. Preserve source activity/assay/target/document IDs.
No automatic cross-study averaging, endpoint conversion or categorical
selectivity labels. Ratios require compatible positive exact measurements on
an appropriate ratio scale. Record arithmetic direction and interpretation.

## Context QC has two separate scopes

1. `pharmacology_context_match_status`: match **between the two drugs' assays**.
   `exact_context_match` means the recorded comparison key agrees, not that
   unreported experimental conditions were independently verified.
2. `structure_pharmacology_context_match_status`: match **between the pharmacology
   receptor and the structural block**. Never substitute the assay-to-assay
   status for this status.

`exact_context_match`: relevant recorded identity/context fields agree and the
structural construct is verified where applicable.
`partial_context_match`: no explicit conflict but missing conditions, unresolved
splice form or a provisional local construct prevent exact matching.
`context_mismatch`: a known species/composition/context difference exists.
A known mismatch takes precedence over partial/unknown evidence.

**β2 pharmacology versus β3 IFP is context_mismatch**, even if the two drugs'
pharmacology shares an assay and document. The 9CTJ-derived α2 C/D/E construct
is provisional, not a matched full α2β3γ2 pentamer. Preserve these constraints.

## Interpretation and immutable provenance

Two drugs permit descriptive contrasts, not evidence that a feature determines
affinity/efficacy or predicts selectivity, phenotype or adverse effects.
Separate supported observations, untested hypotheses and unsupported claims.
Drug poses are not independent experimental replicates. Do not treat repeated
integration rows as additional observations. Avoid post hoc causal scoring.

Read and cite this contract before each new analysis. Snapshot its version/hash
in the run. Preserve all previous runs and inputs, including LATEST_RUN.txt when
existing-file immutability is requested. Save new results in a new run with code,
input provenance, tests and QC. A changed scientific axis requires an explicit
contract revision rather than silently changing the primary analysis.
