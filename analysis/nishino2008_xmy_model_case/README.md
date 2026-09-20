# Nishino 2008 X–M–Y model case

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
