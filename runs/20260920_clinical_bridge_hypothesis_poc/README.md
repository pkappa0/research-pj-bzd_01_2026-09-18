# Clinical bridge hypothesis PoC — 2026-09-20

Read reports/FINAL_REPORT.md first. This run extends the analysis contract as explicitly requested; source contract and user request are archived under config/.

Stages: `python code/acquire.py` freezes unchanged X, writes prespecified clinical rules, retrieves cached raw API responses, builds/QCs Y and writes CLINICAL_RESULTS_FROZEN.json. `python code/analyze.py` verifies freeze hashes before correspondence and plots. `python code/validate.py` independently checks contingency arithmetic, missingness, freeze hashes, provenance and dimensions. Run from any directory with Python plus numpy/pandas/scipy/matplotlib/requests. No API key required for saved query forms. Raw responses permit cached reproduction; to query a new release create a new run, never overwrite this run's data.

Primary estimand is all-role openFDA report-document reporting disproportionality, not primary-suspect disproportionality, incidence, drug effect or risk. The API does not expose FAERS PS/SS distinctions. A report-level AND with drugcharacterization would incorrectly include roles belonging to co-reported drugs, so no such sensitivity is presented. Exact harmonized generic OR exact medicinalproduct is used for every drug; unknown brands/salts are not force-mapped. Co-reported brands in drug-array marginal counts are not drug-identity synonyms. ChEMBL synonym lists are separately sourced and are not added into report totals.

Demographic marginals retain sex, reported age, country, receivedate and seriousness; top-100 distributions may truncate (especially dates). Age values are not converted to years because age units are not jointly available. These are not patient-level or event-specific covariates. Query response examples are raw API examples, not a complete microdata cohort. No demographic adjustment, independent case deduplication, or role-specific analysis is claimed.

Plots retain original feature ordering. Alpha1 and alpha2 blocks remain identifiable; no structural clustering or feature selection was performed. A maximum absolute correlation is descriptive across all 672 comparisons and is not a validated discovery.
