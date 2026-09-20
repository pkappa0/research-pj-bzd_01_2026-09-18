from pathlib import Path
import json,hashlib,shutil,subprocess,sys,importlib.metadata as md
import pandas as pd
R=Path(__file__).resolve().parents[1];ROOT=R.parents[1]
summary=json.loads((R/'qc/model_case_summary.json').read_text());commitfile=R/'config/results_commit.txt';commit=commitfile.read_text().strip() if commitfile.exists() else 'RESULTS_COMMIT_TO_BE_STAMPED'
report=f'''# 1. Research question

Can Structural X, controlled in vivo M and Clinical Y be aligned for the same four drugs in a reproducible representation? **Yes, with explicitly partial clinical coverage and an incompletely specified reported ED50 calculation context.** X is complete for four drugs; all four reported M values are available; Y contains 14/16 eligible cells, with three PTs complete across four drugs. This is a representation-establishment model case, not a predictive or causal validation.

# 2. Why Nishino 2008 was selected

[Nishino et al. 2008, J Pharmacol Sci 107:349–354, DOI 10.1254/jphs.08107FP](https://www.jstage.jst.go.jp/article/jphs/107/3/107_08107FP/_pdf) provides multiple drugs in one controlled rotarod protocol, a dose/time count table and reported probit ED50 values. The paper itself was checked, including visual review of Table 1. This creates a stronger common in vivo context than pooling unrelated assay summaries, while preserving the paper's reporting limitations.

# 3. Drug-set definition

Primary fixed order: diazepam, triazolam, brotizolam, lormetazepam. Chemical identities are linked to ChEMBL IDs and source SMILES. Rilmazafone remains in the complete five-drug raw M table and ED50 table but is excluded from primary structural/matched analyses because metabolic activation produces active species. See docs/rilmazafone_exclusion_note.md. Parent docking must not be equated to activity after oral prodrug administration. Direct-acting parent status does not eliminate metabolism or exposure confounding for the primary four.

# 4. Structural X

X remains four independent views: PLIF_alpha1, PLIF_alpha2, Vina_alpha1 and Vina_alpha2. No composite or PLIF weighting is created.

Environment audit verified Vina binary, scientific package versions, conda builds, original code and PLIP effective settings. Replaying existing ligand preparation reproduced diazepam/triazolam ligand bytes. Receptor scientific records reproduced; only path-bearing headers differed. Original frozen receptor bytes were used in production. The pip installer version difference has no effect on the verified calculations. Therefore existing frozen diazepam/triazolam outputs were reused and only brotizolam/lormetazepam were added. Details, hashes, boxes and preparation rules are in qc/nishino_structural_environment_check.md.

All 40 expected drug/receptor/seed outputs are present. All 354 returned poses are retained (175 reused, 179 new). New PLIP extraction succeeds for all 179 poses. Every output PDBQT matches its docking-manifest hash; all rank-1 poses agree with seed minima. All eight receptor summaries have five seeds. Primary score is the median seed minimum; mean, sample SD, IQR, min/max and negative-score favorability are stored separately.

| Drug | alpha1 median Vina docking score | alpha2 median Vina docking score |
|---|---:|---:|
| diazepam | -10.130 | -7.774 |
| triazolam | -10.887 | -6.820 |
| brotizolam | -10.215 | -6.968 |
| lormetazepam | -8.857 | -6.876 |

Units are kcal/mol as emitted by Vina, not experimental affinity or true binding free energy. Receptor raw scales are not interpreted as equivalent. Lormetazepam/alpha1 has a seed score of -8.077 versus four near -8.86, retained without exclusion; sample SD is 0.349. Seed repeats are computational, not biological replicates.

All 48 frozen PLIF features remain, in original order: 21 alpha1 and 27 alpha2. Frequency uses the union of contact-positive poses across interaction types divided by evaluable poses. New contacts outside the frozen schema are flagged, not added: one brotizolam/alpha1 contact at BZD_SITE_159 SER, and two lormetazepam/alpha2 contacts at BZD_GAMMA2_191 GLY. These are three interaction records, not three drugs/features. Unmapped contacts: zero. The retained representation is technically compatible but not exhaustive for all new contacts. Twelve of the 48 features are constant across this four-drug subset; they remain in X and receive NA correlations.

Alpha2 remains a provisional 9CTJ local construct. Structural species/composition/chain context is not asserted to match the mouse assay; beta2 pharmacology versus beta3 IFP remains a context mismatch. Lormetazepam is one deterministic configuration from unspecified source stereochemistry. The extra hydroxyl torsion versus RDKit heavy-atom rotatable-bond count was reviewed and documented rather than manually altered.

# 5. in vivo M

Source verification covers protocol, all 80 five-drug dose/time cells plus four separate controls, and all five reported ED50s. The primary four contribute 64 cells. Male ICR mice received oral drug in 0.5% carboxymethylcellulose sodium; 15 rpm, 180-second trials, positive event = fall. Successful completion of two consecutive pretests 24 hours earlier was required; measurement times are 15/30/60/90 min; each cell has ten tested mice. Independent total animal count is not inferred from repeated cells.

| Drug | Reported ED50 mg/kg | Reported interval | −log10(ED50) |
|---|---:|---|---:|
| diazepam | 3.11 | 2.53–3.81 | -0.493 |
| triazolam | 1.25 | 1.00–1.49 | -0.097 |
| brotizolam | 5.76 | 5.45–9.12 | -0.760 |
| lormetazepam | 3.39 | 2.75–4.27 | -0.530 |

M_raw is reported ED50. The logarithmic sign reversal is visualization only. ED50 spans 4.608-fold, but precision and protocol limits matter. The paper says probit, without specifying the ED50 time aggregation or interval confidence level. These are reported intervals, not asserted 95% CIs. No primary refit is performed. The printed triazolam 0.5 mg/kg/90 min value is 9/10; it remains flagged and unchanged. Rotarod motor-impairment potency is not clinical risk, incidence or binding affinity. See docs/nishino_source_and_protocol_notes.md.

# 6. Clinical Y

The existing openFDA/FAERS report-document method is reused: receivedate 2004-01-01 through 2025-12-31; all drug roles; exact uppercase generic-name OR exact medicinalproduct name; no brand expansion or alias-count summation. C1 grouping is provisional; the actual variables are the four PT strings. No new term selection, alternative drug-specific database or missing-value imputation is introduced.

Background total 20,018,532 and all four PT background counts match the frozen snapshot. Fresh counts for diazepam/triazolam also match, permitting frozen Y reuse. New compounds use the same four-cell ROR formula, 0.5 correction to all cells only when a zero exists, and the same a<5 primary mask. API responses are saved as counts/metadata, with exact query URL, timestamp and hash; incidental individual reports are not archived. No patient-level deduplication or covariate adjustment is claimed.

| Drug | Total target reports | Ataxia a | Balance disorder a | Coordination abnormal a | Gait disturbance a | Available Y |
|---|---:|---:|---:|---:|---:|---|
| diazepam | 118049 | 375 | 1355 | 335 | 2287 | 4/4 |
| triazolam | 3893 | 12 | 29 | 4 | 62 | 3/4 |
| brotizolam | 4473 | 7 | 6 | 0 | 38 | 3/4 |
| lormetazepam | 5825 | 13 | 65 | 8 | 48 | 4/4 |

Raw estimability and primary eligibility are distinct fields in the audit. Triazolam/brotizolam Coordination abnormal are NA after the original low-count rule. Brotizolam Ataxia and Balance disorder pass at only 7 and 6 reports; passing the mask is not a guarantee of robust clinical evidence. Frozen/raw a,b,c,d and comparator totals are preserved. logROR describes reporting disproportionality, neither incidence nor a risk ratio. Negative logROR values remain visible on the diverging clinical heatmap.

# 7. Matched X–M–Y representation

`data/nishino_primary4_XMY_master.csv` has one row per drug with receptor-specific PLIF references, raw Vina medians and favorability, reported ED50/interval/potency, four PT logROR values and layer completeness. All joins use canonical drug identity. Structural sources, assay values and clinical sources are not numerically collapsed. No causal arrows appear in the synthesis figure. Complete M denotes reported-value availability; method uncertainty remains explicit.

# 8. Descriptive correspondence

The observational unit is drug, n=4 at most. No p-values, significance-based selection, residue ranking or learned model are used. For this user-requested n=4 case only, a documented floor of three complete drugs is used; constant inputs and n<3 give NA. Historical ten-drug correlation tables and n>=5 rules remain unchanged.

X ↔ M: Vina favorability versus rotarod potency has rho +0.4 for alpha1 and -0.4 for alpha2, each n=4. All 48 PLIF feature correlations are stored in frozen order; 36 are estimable and 12 constant inputs are NA. Opposite signs are reported without receptor-dominance or superiority claims.

M ↔ Y: Ataxia rho +0.8, Balance disorder +0.4, Gait disturbance +0.6 (each n=4). Coordination abnormal has n=2 and rho NA. These values describe ordering in this selected four-drug case; they do not validate causality or predict patient outcomes.

X ↔ Y: 200 secondary rows (48 PLIF + 2 Vina views × four PTs) are retained with n/status, including NA for insufficient coverage or constant inputs. They do not replace the larger ten-drug PoC. No positive/negative result is used to alter X, choose features or define success.

# 9. What this model case establishes

A reproducible identity-aligned multi-view X, common-study reported M and explicitly masked Y can be assembled for these four drugs. Existing structural outputs can be extended after a meaningful environment replay rather than assumed compatible. Traceable raw values, out-of-schema contacts, low-count clinical cells and source anomalies remain accessible.

# 10. What it does NOT establish

No mediation, causal residue effect, receptor dominance, clinical prediction, experimental affinity, general representation superiority or model validation is established. Four drugs do not support supervised ML. A common study does not remove pharmacokinetic confounding or unreported assay-summary details. Repeated poses, seeds and time points do not increase independent drug n.

# 11. Main bottleneck identified

The limiting layers are endpoint-specific clinical coverage and in vivo summary-method documentation, not basic structural computation or identity linkage. The four-PT Y is incomplete, although three PTs have complete four-drug coverage. Some eligible cells are still sparse. Nishino's reported ED50 method needs clarification before strong cross-study numeric harmonization. Frozen structural representation additionally omits the explicitly flagged new contacts and samples one unspecified lormetazepam configuration.

# 12. Go / no-go reasoning

See docs/nishino_model_case_go_nogo.md. CASE 1 is supported with qualifications for representation and systematic matched-study expansion; CASE 2 applies locally to Coordination abnormal coverage. CASE 3 is not an overall structural harmonization failure, but M method details should be resolved before pooling potency estimates. None of these decisions depends on correlations being positive.

# 13. Next experimental / data milestone

First clarify ED50 time basis and uncertainty definition and expand Tier-1 same-study multi-drug in vivo evidence. Audit clinical PT coverage using a consistent source/protocol before selecting a larger drug panel. Preserve exact dose/time grids as context-specific observations. A new in vivo experiment is not automatically justified by this model case; future prioritization needs validated uncertainty and explicit assay design.

# 14. Future ML role

Only a roadmap is provided: schema → matched in vivo collection → expanded drug n → possible receptor expansion → prespecified Vina-only, PLIF-only, combined X and X+available-M predictive baselines → sparse-M methods → external drug validation → possible uncertainty/information-gain experiment prioritization. No algorithm is selected or trained here. See docs/future_ml_roadmap.md and docs/future_in_vivo_expansion_plan.md.

# 15. Reproducibility

Analysis directory: `analysis/nishino2008_xmy_model_case/`. Source hashes/configurations, exact source requests, prepared ligand/receptor files, four compact per-drug archives, raw tables, independent tests and figure scripts are retained. Original 9/20 runs and LATEST_RUN.txt are unchanged. The complete artifact index is OUTPUT_INDEX.md; reproduction commands are in README.md.

Six independent validation tests pass: all archive score/RMSD values and seed summaries; new PLIF frequencies rebuilt from archived XML; existing-drug equality/frozen hashes; clinical counts/ROR mask; independent rank-Pearson check of X–M; and raw M/identity alignment. All four figure sets were visually reviewed. qc/validation_log.txt records the run. No historical file was modified.

Final analysis/results Git commit: `{commit}`. A subsequent documentation-only commit records this immutable results hash, avoiding a self-referential Git hash inside its own commit. Resolve the documentation tip with `git log -1 --format=%H -- analysis/nishino2008_xmy_model_case/FINAL_REPORT_NISHINO_XMY_MODEL_CASE.md`.
'''
(R/'FINAL_REPORT_NISHINO_XMY_MODEL_CASE.md').write_text(report)
readme='''# Nishino 2008 X–M–Y model case

Start with FINAL_REPORT_NISHINO_XMY_MODEL_CASE.md and OUTPUT_INDEX.md. Primary result: four-drug X–M is populated; Clinical C1 is available in 14/16 cells, with three PTs complete across all four. Reported ED50 context remains partly unspecified. No ML, composite predictor or causal model.

## Reproduce

Run from repository root. Structural stages require the exact original layered environment, recorded in qc/environment_check.json. Python executable on this machine:

`/Users/k-atsumi/Documents/Codex/2026-09-18/github-readme-final-report-csv-git/work/standardized_env/venv/bin/python`

Use a writable MPLCONFIGDIR. Plot/clinical/report stages require numpy, pandas, scipy, matplotlib and requests. Core scripts import the unchanged repository src modules, whose byte-identical source snapshots are provided. Full reproduction requires the repository and original standardized source run, not just this exported folder. The original paper remains at its publisher URL and local verification path in QC.

Order (substitute the exact environment executable for python for structural stages):

```sh
python analysis/nishino2008_xmy_model_case/scripts/acquire_identity.py
python analysis/nishino2008_xmy_model_case/scripts/check_environment.py
python analysis/nishino2008_xmy_model_case/scripts/structural_pipeline.py prepare
python analysis/nishino2008_xmy_model_case/scripts/structural_pipeline.py dock
python analysis/nishino2008_xmy_model_case/scripts/structural_pipeline.py plip
python analysis/nishino2008_xmy_model_case/scripts/assemble_structural.py
python analysis/nishino2008_xmy_model_case/scripts/extract_in_vivo.py
python analysis/nishino2008_xmy_model_case/scripts/acquire_clinical.py
python analysis/nishino2008_xmy_model_case/scripts/analyze_model_case.py
python analysis/nishino2008_xmy_model_case/scripts/plot_model_case.py
python analysis/nishino2008_xmy_model_case/scripts/validate_model_case.py
python analysis/nishino2008_xmy_model_case/scripts/write_documentation.py
python analysis/nishino2008_xmy_model_case/scripts/finalize_report.py
```

The structural code has no M/Y input. Existing frozen diazepam/triazolam poses and PLIF are reused; only new compounds are computed. Docking/PLIP checkpoints reside outside the repository at work/nishino_jobs. Do not delete caches and redock merely to match outcomes. New outputs are archived by drug with PDBQT, pose SDF, XML/TXT, logs and metadata; large reconstructable complexes are omitted using the existing convention. Old drugs' archives are copied byte-for-byte. Clinical HTTP errors are not counted as zero; only explicit NOT_FOUND responses imply zero API matches. The script fails rather than mix a changed background with frozen Y.

Feature definitions and correlation floor are in config; source/code/manifests and QC retain deviations. Minor runtime timestamps in regenerated provenance may change, while structural values and scientific source snapshots should remain stable. The current task is authorized to compute n=4 descriptive correlations and new-drug structural data; historical no-recompute visualization tasks remain untouched.
'''
(R/'README.md').write_text(readme)
for src,dest in [('nishino_environment.log','environment_replay.log'),('nishino_prepare.log','ligand_receptor_prepare.log'),('nishino_dock.log','new_docking.log'),('nishino_plip.log','new_plip.log'),('nishino_clinical.log','clinical_acquisition.log')]:
 p=ROOT.parent/src
 if p.exists():shutil.copy2(p,R/'qc'/dest)
(R/'config/analysis_software_versions.json').write_text(json.dumps({'python':sys.version,**{p:md.version(p) for p in ['numpy','pandas','scipy','matplotlib','requests']}},indent=2)+'\n')
# Explicit source linkage to original independent layers and unchanged mapping.
sources=['ANALYSIS_CONTRACT.md','runs/20260919_0724_10drug_standardized_redocking/config/docking_config.json','runs/20260919_0724_10drug_standardized_redocking/tables/docking_run_manifest.csv','runs/20260919_0724_10drug_standardized_redocking/tables/fingerprint_10drug_multireceptor.csv','runs/20260917_1706_structure_mouse_bridge/results/tables/corrected_residue_mapping.csv','runs/20260920_clinical_bridge_hypothesis_poc/tables/clinical_fingerprint_10drug.csv','runs/20260920_clinical_bridge_hypothesis_poc/tables/clinical_disproportionality.csv','runs/20260920_clinical_bridge_hypothesis_poc/CLINICAL_RESULTS_FROZEN.json']
(R/'config/upstream_source_manifest.json').write_text(json.dumps([{'path':p,'sha256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in sources],indent=2)+'\n')
purposes={'data':'Reproducible source, derived or matched data table','qc':'Quality-control evidence or verification log','config':'Fixed method, software version or provenance record','scripts':'Reproduction or verification code','docs':'Scientific interpretation, source limitations or future plan','tables':'Structural preparation, mapping or extraction audit','figures':'Publication figure','raw':'Raw chemical response, prepared structural input or archived computational evidence'}
lines=['# Output index','','Paths are relative to this analysis directory. All generated outputs are listed; no individual animal/patient inference is made.','']
for p in sorted(R.rglob('*')):
 if not p.is_file() or '__pycache__' in p.parts or p.name=='OUTPUT_INDEX.md':continue
 rel=p.relative_to(R);purpose=purposes.get(rel.parts[0],'Entry-point report or reproduction guide')
 if rel.parts[0]=='figures':purpose={'nishino_primary4_XMY_model_case':'Identity-aligned PLIF, Vina, reported M and masked Y','nishino_primary4_rotarod_profiles':'Raw dose/time profiles and separate reported ED50 intervals','nishino_primary4_X_M_correspondence':'Frozen-order descriptive structural correspondence with M, n=4','nishino_three_layer_summary':'Noncausal concept and empirical layer completeness'}[p.stem]
 lines.append(f'- [{rel}]({rel}) — {purpose}.')
lines.append('- [OUTPUT_INDEX.md](OUTPUT_INDEX.md) — Complete generated-output inventory.')
(R/'OUTPUT_INDEX.md').write_text('\n'.join(lines)+'\n')
print('Report, README, provenance, and complete output index written.')
