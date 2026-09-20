from pathlib import Path
import json,hashlib,subprocess,importlib.metadata as md,sys
import pandas as pd
R=Path(__file__).resolve().parents[1];ROOT=R.parents[1]
def write(p,s):(R/p).write_text(s.strip()+'\n')
write('config/CONTRACT_ADDENDUM.md','''# Authorized Nishino matched-model-case extension

The user's MASTER TASK authorizes an isolated four-drug representation study, including new standardized docking for absent compounds, unchanged 48-feature PLIF extraction, the frozen clinical ROR methodology and descriptive n=4 X–M–Y correspondence. This supersedes earlier visualization-only/no-redocking scope only for the present new run. Historical results, original analysis contract, source feature selection and receptor blocks are unchanged. No supervised ML, composite feature, mediation or outcome-based optimization is authorized or performed. The main comparison remains drugs within a fixed receptor. Refer to docs/USER_REQUEST.md for the complete task.

The n=4 descriptive correlation task uses a documented minimum of three complete drugs; two-drug correlations are suppressed as NA. This run-specific descriptive floor does not modify the existing ten-drug n>=5 analysis. Clinical count masking remains the original a<5 rule.
''')
write('docs/rilmazafone_exclusion_note.md','''# Rilmazafone remains in the source data, outside the primary model case

Nishino 2008 includes rilmazafone. All 16 dose/time cells and its reported ED50 9.55 mg/kg with interval 8.20–12.9 are retained in the five-drug raw tables. Rilmazafone is excluded from the primary four-drug structural, matched-master and correspondence tables.

It is a prodrug that generates pharmacologically active cyclic metabolites. Docking the administered parent cannot represent the oral pharmacological activity after metabolic activation. This is a compound-identity/exposure issue, not grounds to discard the source results. A later metabolism-aware sensitivity case would need metabolite identities, exposure/time data and separate parent/metabolite structural representations; none is substituted here.

Sources: [Nishino 2008](https://doi.org/10.1254/jphs.08107FP), [Koike et al., metabolite structure determination, 1988](https://pubmed.ncbi.nlm.nih.gov/3381539/), [active cyclic metabolite transport, 1993](https://pubmed.ncbi.nlm.nih.gov/8499579/). The four primary compounds are treated as direct-acting parents, but this does not imply that their in vivo effects are free of pharmacokinetic or metabolite contributions.
''')
env=json.loads((R/'qc/environment_check.json').read_text())
write('qc/nishino_structural_environment_check.md','''# Structural environment check and decision

**Decision: reuse frozen diazepam/triazolam; add brotizolam/lormetazepam only.** No material computational-environment uncertainty remained after the checks below. No incompatible batch is mixed.

- Exact Vina 1.2.7 executable SHA-256: `823c2bbacf26d72183861322345f0a89736aca66c8e81054c66f93af5ad623f1`.
- Python 3.12.13, RDKit 2025.9.6, Meeko 0.7.1, PLIP 3.0.0; original layered environment retained at the executable path recorded in environment_check.json. All original conda package versions/builds/hashes matched. Scientific Python distributions matched. Only pip differs (26.1.1 recorded vs 25.0.1 current); this installer is not used to calculate structures and is immaterial to the verified replay. No environment packages were altered.
- Original preparation code, extraction code and effective PLIP config match the frozen source bytes. PLIP NOHYDRO/NOFIX/NOFIXFILE are true and interaction thresholds are unchanged. Read source_standardized_plip.py for atom reconstruction and XML extraction.
- Replayed diazepam and triazolam ligand SDF/PDBQT are byte-identical. Receptor ATOM/HETATM coordinates, atom types, charges, connectivity and all scientific records reproduce exactly. Only path-bearing COMPND/REMARK Name headers differ. Final production uses the original frozen receptor files verbatim, with original hashes, not replay output headers.
- Alpha1 receptor PDBQT SHA-256: `1089f81e8bbc900a2522d34c5ed65a11c35f3c8046954386dafde623056a806a`.
- Alpha2 receptor PDBQT SHA-256: `c3bcd380238909d2b6c7cabc56ab6cb69e4f80613e7640cf760aa1f4a737ab87`.
- Boxes are fixed **within receptor** across drugs; alpha1 and alpha2 do not have identical coordinate centers. Alpha1 center (117.44499969482422,157.46474990844726,110.48740043640137), alpha2 (143.42579952820148,96.86254491988079,139.9939981778101); both 20×20×20 Å. Alpha2 is the original provisional transferred box, not a new independently optimized pocket.
- Exhaustiveness 8; num_modes 9; energy_range 4; cpu 1 per job; min_rmsd 1; scoring vina; seeds 2026091901–2026091905. All returned poses retained, no padding or score-based reruns.
- Ligands use source ChEMBL SMILES/formal charge/tautomer, explicit H, ETKDGv3 seed 20260919, one thread, enforce chirality, no random coordinates, MMFF94s max 2000 iterations, convergence required; Meeko Gasteiger with rigid macrocycles and nonflexible amides. No pH/tautomer enumeration. Lormetazepam source stereochemistry is unspecified: a deterministic sampled configuration is not a racemate ensemble or experimentally assigned enantiomer.
- Receptors originate from archived ATOM-only selected chains, model 1, altloc blank/A, no water/heteroatoms, OpenBabel add-H and Gasteiger without pH titration. OpenBabel emits an aromatic-kekulization warning during replay; it is retained in the log and identical scientific-record replay was verified. Frozen prepared receptor bytes are the production definition.
- Residue map is the existing corrected map. All frozen F01–F48 identities remain valid. Any-contact frequency is the union of contact-positive pose IDs across interaction types divided by all successfully extracted returned poses, pooled across seeds. No sum of marginal interaction frequencies.

See environment_check.json, source replay logs, config/standardized_* locks, tables/residue_mapping_qc.csv, config/PREPARED_INPUTS_FROZEN.json and qc/docking_file_manifest.csv. Alpha2 9CTJ C/D/E remains a provisional local construct; β2 pharmacology versus β3 IFP is a context mismatch. Mouse assay context is not asserted to match these structural preparations.
''')
write('docs/nishino_source_and_protocol_notes.md','''# Direct source verification and M limitations

The authoritative source is [Nishino et al. 2008, J Pharmacol Sci 107:349–354](https://www.jstage.jst.go.jp/article/jphs/107/3/107_08107FP/_pdf), DOI 10.1254/jphs.08107FP. The locally retained original PDF was read and Table 1 rendered for visual checking; its exact hash/path are recorded in qc/nishino_source_verification.json. Values were checked against the paper, not inherited from chat text. A licensed full-paper redistribution is not needed to reproduce the extracted numeric tables; the publisher link identifies the source.

Male ICR mice, 23–28 g, received oral drug in 0.5% carboxymethylcellulose sodium (10 mL/kg). The rod rotated at 15 rpm and was 3 cm in diameter. Selection used two successful consecutive trials 24 hours before testing. Failure meant falling within the 180-second trial. Observations were at 15, 30, 60 and 90 minutes. Each table cell has ten tested mice. The total number of independent animals and repeated-measure allocation are not explicit; 80 cells are not 800 independent animals. The five-week age statement belongs to the plus-maze procedure and is not silently assigned to rotarod.

Five drugs × four doses × four times = 80 raw response cells; four vehicle-control cells (0/10) are stored separately. The primary four contribute 64 cells. Triazolam 0.5 mg/kg at 90 min is printed as 9/10 without a significance marker. This unusual entry is retained and flagged, not corrected to 0/10 or another inferred value.

Reported ED50 is primary; the paper states probit estimation. Results list diazepam 3.11 [2.53,3.81], triazolam 1.25 [1.00,1.49], brotizolam 5.76 [5.45,9.12], lormetazepam 3.39 [2.75,4.27] mg/kg; rilmazafone 9.55 [8.20,12.9]. The task's ed50_ci_low/high fields store those parenthetical bounds. The source does not explicitly state their confidence level or precisely which time/aggregation enters the probit estimate. They are therefore labeled reported intervals, not asserted 95% CIs. The Results mention diazepam 5/10 at 60 minutes after 5 mg/kg immediately before ED50=3.11; no attempt is made to repair or refit this summary from Table 1. Clarifying the ED50 basis is a useful next source-verification milestone.

M_raw remains ED50 in mg/kg; rotarod_potency = −log10(ED50_mg_kg) is only a monotonic visual representation. It is neither clinical incidence/risk nor binding affinity. Lower ED50 indicates higher potency under this assay. Systemic exposure, metabolism, sedation and motor effects can contribute; receptor selectivity or mediation cannot be inferred. No ED50 re-estimation is used in this task.
''')
write('docs/nishino_model_case_go_nogo.md','''# Go / no-go reasoning: a qualified representation success

The conclusion is driven by coverage, traceability and harmonization, not correlation sign.

| Process criterion | Observed result | Implication |
|---|---|---|
| A Structural feasibility | 4×48 frozen PLIF; 8 receptor summaries; 40 seeds; 354 poses | Feasible under a verified common environment; retain construct and single-conformer caveats |
| B In vivo feasibility | Four reported ED50s and intervals; 64 primary raw response cells | Common-study M is available, but ED50 time basis and interval level need clarification |
| C Clinical feasibility | 14/16 usable PT cells; three PTs complete across all four drugs | Broad sparse-drug failure did not occur; Coordination abnormal remains a targeted coverage gap |
| D Alignment feasibility | Unique drug identity joins all four rows without imputation | Matched representation is technically constructible |
| E Dynamic range | 36/48 PLIF features vary; 12 are constant; ED50 1.25–5.76 mg/kg (4.608-fold); receptor medians and available Y vary | Sufficient observed variation to motivate collecting more matched data, not sufficient to establish predictive information |
| F No circularity | Fixed geometry, seeds, schema and preparation; no M/Y-driven selection or score tuning | X is computationally independent of M/Y construction |

**CASE 1, qualified:** proceed to systematic matched in vivo data collection and protocol audit. Three prespecified C1 PTs are complete for all four drugs; the representation is useful despite incomplete fourth-PT coverage. This is not permission to deploy or train an n=4 predictive model.

**CASE 2, localized:** Coordination abnormal has target-event counts 4 for triazolam and 0 for brotizolam, below the existing five-report threshold. The gap is clinical coverage, not evidence of absent biological effects. An appropriate clinical-source validation may later be needed, using a common protocol for all drugs; no individual-drug database substitution is made here.

**CASE 3, not an overall structural failure:** the fixed feature matrix is harmonizable. Three mapped contact records fall outside the 48-feature source schema: one brotizolam/alpha1 BZD_SITE_159 SER and two lormetazepam/alpha2 BZD_GAMMA2_191 GLY; they are explicitly preserved but do not expand it. M's reported summary method remains incompletely described. Resolve that assay-summary uncertainty before pooling ED50 across studies or making stronger quantitative M claims.

Next milestones: clarify Nishino ED50 calculation context where possible; retain the exact dose/time endpoint as fallback context; gather additional Tier-1 within-study multi-drug assays; audit PT-specific clinical coverage across candidate drugs before expanding X. Receptor expansion and ML remain later-stage work.
''')
write('docs/future_in_vivo_expansion_plan.md','''# Future in vivo expansion plan

Keep one raw observation per drug × study × assay context × dose × time × endpoint, linked to an immutable source location. Required metadata: canonical parent/active-moiety identity and identifiers; species; strain; sex; age/weight where available; route; formulation/vehicle; dose and unit; observation time; assay protocol; rotarod speed/acceleration and rod geometry; training/selection; positive-event or latency definition; n and repeated-measure design; effect value; uncertainty/CI type and level; reported summary-model method; source page/table and extraction QC.

Explicitly distinguish administered prodrug from active metabolites. Preserve reported units and raw counts; do not pool latency with binary failure or equate different dose/time windows. Record absent metadata as unavailable. Include source anomaly flags and double-check nonmonotonic values without silently correcting them.

- Tier 1: multiple drugs in one study under the same experimental context. Priority acquisition tier. Keep study as a context block and compare drugs within it.
- Tier 2: closely matched protocols across studies with compatible species, route, speed, training, endpoint and time. Require a documented comparability audit and heterogeneity assessment before any combined model.
- Tier 3: heterogeneous protocols or insufficient metadata. Retain as qualitative/contextual evidence; do not pool numerical values as equivalent observations.

Nishino demonstrates why reported ED50 alone is insufficient metadata. Collect the dose/time grid, reported probit inputs/time basis and interval definition alongside it. The primary selection criterion for new studies is context completeness and drug overlap, not favorable or positive correlations. No additional animal experiments are performed or automatically proposed as validated designs here.
''')
write('docs/future_ml_roadmap.md','''# Future ML roadmap — planning only

No ML is implemented in this model case.

0. Current: establish X/M/Y schema, source provenance, independent structural views and explicit missingness.
1. Systematically collect matched in vivo datasets, prioritizing Tier-1 multi-drug studies.
2. Expand independent drug n; retain raw Clinical PT-level outcomes, provisional C1 grouping and common coverage rules. Do not replace missing PTs with a score.
3. Expand receptor representation (e.g. alpha3/alpha5) only after sufficient independent drugs and reproducible context exist; preserve receptor identity.
4. Compare prespecified predictive baselines: A Vina-only → Y; B PLIF-only → Y; C PLIF+Vina → Y; D PLIF+Vina+available M → Y. Decide outcomes, drug-level held-out evaluation, simple baseline, uncertainty and leakage controls before training. This future plan does not imply any current superiority or causal chain. Keep all poses/seeds of a drug in one split; fit preprocessing only on training drugs. Account for scaffold/study dependencies and reserve external drugs.
5. Investigate sparse-M approaches: multi-task, multi-view, missing-label/partially observed intermediate-phenotype methods and latent representations. Do not select algorithms before coverage and validation requirements are known; explicitly evaluate missingness assumptions.
6. External drug validation under a frozen protocol and endpoint definitions; report failures and calibration, not only successful correlations.
7. Only if validated, explore uncertainty/information-gain-based prioritization of new in vivo characterization. Active-learning or experimental-design criteria must be evaluated against cost, context compatibility and prospective validation.

n=4 is unsuitable for supervised model selection. Clinical reporting signals remain observational logROR, not incidence or risk ratios. An intermediate measured phenotype does not establish mediation. Sample size sufficiency must be evaluated for the eventual target/model and independent validation plan, not inferred from the 48 feature count.
''')
