# BZD-site pharmacology triage

`src/bzd_pharmacology_filter.py` classifies the 434 activities already acquired in
`runs/20260918_1117_chembl_acquisition`. It makes no network requests, imputes no
values, and does not alter existing runs, IFP/phenotype inputs, or LATEST_RUN.txt.
Configuration: `config/bzd_pharmacology_filter.json`.

## Run

```sh
MPLBACKEND=Agg .venv/bin/python -m pytest -q --junitxml=/tmp/bzd-tests.xml
.venv/bin/python -m src.bzd_pharmacology_filter \
  --config config/bzd_pharmacology_filter.json --test-results /tmp/bzd-tests.xml
```

Every invocation creates a new immutable `runs/*_bzd_pharmacology_filter/`.
Classification is deterministic; timestamps, pre-run Git commit and test timings
are provenance rather than analytical inputs. A failed QC stops the run and,
if its directory has already been created, marks its manifest FAILED.
The provided JUnit XML records a prior test invocation; run it against the current
code before generating an official run. Unit/integration tests run offline.

## Interpretation

- PRIMARY requires explicitly matching full single-alpha1/alpha2, single-beta,
  gamma2 composition in both target name and assay text; specific PROTEIN COMPLEX;
  meaningful endpoint; no review flags. This conservative policy is stricter than
  accepting a subtype inferred from one annotation alone.
- SECONDARY retains generic, single-subunit, incomplete, native or other-subtype
  context. Group members in target metadata are not treated as coassembled subunits.
- REVIEW takes precedence over every other category. Subunit, organism, B/F
  annotation or endpoint conflicts are not repaired, even when host-cell species
  or coarse assay annotations plausibly explain them. Original papers are needed.
- EXCLUDE retains out-of-scope rows and an explicit reason, including non-BZD
  probe/channel readouts or whole-animal doses. No scaffold/name-based exclusion.

The classifier is a conservative rule-based triage, not a validated language model
or exhaustive assay parser. Unhandled syntax, negation and complex constructs
remain limitations. A PRIMARY label does not establish direct BZD-site binding or
experimental comparability by itself. Receptor-composition support for BZD relevance
is explicitly labelled `explicit_abg_receptor_candidate`, distinct from an explicit
BZD-site/probe statement. Original methods and the exact IFP structure still need
checking, including species, beta subunit and gamma splice form.

Target name, target components, assay description and each organism field are
retained separately. Analysis species uses reported assay organism only; missing or
conflicting evidence yields unknown. Description-only species is evidence, not an
imputed species value. Composition is extracted from explicit assay tokens, never
a union with target metadata. Mixed alpha1/alpha2 values are not assigned twice.

Binding IC50 is distinct from functional IC50. Functional efficacy includes response
measurements, with maximal/fixed-concentration/unspecified-response distinctions.
No conversions to common pChEMBL labels, averaging, reference-drug normalization or
cross-endpoint ratios are performed. Numeric strings, nulls, inequalities and source
IDs survive unchanged in the authoritative JSONL; CSV empty cells mean missing and
nested structures are JSON-encoded.

`point_comparison_usable` additionally requires a finite exact value, reported unit,
known assay species and no potential_duplicate annotation. PRIMARY records failing
this condition remain visible. Duplicate flags are not evidence that rows can be
automatically deleted. Reported human data is a coverage indicator, not a reason to
discard rodents or other species.

Comparison groups are compound-independent but assay/document-specific, retaining
species, complete composition, class, endpoint/type/unit/relation and recorded
parameters. Cross-subtype pair candidates additionally require the same compound
and document and an exactly matched recorded context after masking alpha1/alpha2
only. Beta2/beta3 and gamma forms are never masked. Missing conditions are not
established equal; these are candidates for methods checking, not proven matched
experiments. Exact text matching may miss valid pairs. Pair IDs describe the shared
context; the unique pair-table identity is `(molecule_chembl_id, pair_candidate_id)`.

Coverage distinguishes PRIMARY counts from stricter usable counts. Cell-level
manual-review counts use alpha mentions as unresolved hints, never as reassigned
observations; rows with unresolved endpoints are captured in total review counts.
Docking ranks are only assigned among additional compounds with qualifying pair
candidates. Existing IFP drugs are labelled separately. Absence of an eligible new
compound is an acceptable result; missing ranks are not zero scores.

## Outputs and provenance

Each run contains the six requested CSVs and `alpha1_alpha2_pair_candidates.csv`,
a full classified JSONL, FINAL_REPORT.md, QC_REPORT.json, source input snapshots,
code/config/test snapshots and SHA-256 artifact inventory. Original `source_file`
and JSON pointers resolve **relative to the acquisition source run**, not the new
run. `raw/lineage` copies original processed inputs byte-for-byte; original API
response bytes remain in the referenced source run. Original raw responses and
source manifest hashes are verified before classification. All source columns are
retained. `classification_source_run`, `source_processed_line`,
`source_processed_sha256` and `classification_rule_version` identify each lineage.
The manifest Git commit is the HEAD before committing results; snapshot hashes
identify the exact running code, including then-uncommitted additions.

QC checks 1:1 IDs, exact source-field preservation, no REVIEW in PRIMARY, source
ChEMBL IDs, null/inequality retention, no imputed species/composition, strict grouping,
CSV/JSONL roundtrips and unchanged preexisting tracked-file hashes. The test suite
also exercises native/group/single-subunit context, known real annotation conflicts,
endpoint interpretation, beta separation and two runs without mutating prior files.

Policy references: [ChEMBL data FAQ](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions),
[BZD-site structure](https://www.nature.com/articles/s41586-018-0255-3),
[GABAA structural pharmacology](https://www.nature.com/articles/s41586-018-0832-5),
[TBPS channel-site readout](https://pubmed.ncbi.nlm.nih.gov/3035434/).
