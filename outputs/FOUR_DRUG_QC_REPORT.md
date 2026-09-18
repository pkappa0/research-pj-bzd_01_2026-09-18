# FOUR_DRUG_QC_REPORT

QC is descriptive and does not validate biological equivalence of the receptor models.

## Preservation of strict 3-drug outputs

| file | rows |
|---|---:|
| `outputs/docking_results.csv` | 54 |
| `outputs/plip_interactions.csv` | 238 |
| `outputs/fingerprint_binary.csv` | 138 |
| `outputs/fingerprint_frequency.csv` | 138 |
| `outputs/alpha1_alpha2_similarity.csv` | 3 |

## Four-drug computational checks

- Drug × receptor combinations: 8; exactly 9 poses each: **PASS**.
- Vina parameter uniformity (seed=20260917, exhaustiveness=8, num_modes=9, energy_range=4): **PASS**.
- Docking boxes: two recorded boxes from `outputs/docking_boxes.csv`; no result-dependent adjustment: **PASS**.
- Receptor structure IDs: 6HUP, 9CTJ; alpha2 9CTJ construct remains provisional: **PASS**.
- Combined PLIP interaction rows: 339; unresolved mapping rows: 60.
- Fingerprint status counts (including explicit unresolved/no-interaction labels): `{"observed": 83, "unresolved_mapping": 58, "no_interaction_observed": 43}`.

## Biological/QC limitations

- Alpha1 uses 6HUP; alpha2 uses a local 9CTJ-derived construct and a transferred DZP box. They are not an exact matched full-pentamer pair.
- `outputs/residue_mapping.csv` is sequence/structure-alignment based. Mapping absent for a PLIP residue is retained as `unresolved_mapping`; no new common position is inferred.
- `no_interaction_observed` is emitted only for a feature with pose data and resolved mapping. Missing or unresolved evidence is not converted to binary zero.
- Archived ChEMBL rows do not retain assay descriptions for all Tier-1 records, so annotation-vs-assay contradiction checking is explicitly unavailable and no automatic correction is made.
- Alprazolam has only Tier-1 beta2 pharmacology while the fixed docking background is beta3; diazepam includes one beta2 context as well. These are recorded in `results/tables/primary_tier1_qc.csv` as `composition_mismatch_review`.
- No GtoPdb integration, imputation, activity-type/unit mixing, ML or statistical generalization was performed.
