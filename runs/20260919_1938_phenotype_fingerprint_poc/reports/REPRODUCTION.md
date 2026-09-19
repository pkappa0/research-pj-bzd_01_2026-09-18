# Reproduction
Run from repository root using the Python environment recorded by runs/20260919_0724_10drug_standardized_redocking.

The evidence curation is explicit in src/phenotype_fingerprint_poc.py; raw/source_evidence_transcription.csv preserves the transcribed facts. Public URLs, access dates, PDF hashes and extraction locations are in raw/source_manifest.csv.

Set work/phenotype_poc/active_run.txt (relative to this checkout parent) to a NEW run with a copy of config and raw/archived_mouse_phenotype.csv. Then run `python -m src.phenotype_fingerprint_poc`. Existing runs must not be overwritten. No network or docking is used by this analysis step.

Freeze is checked before joining; all 48 columns retain frozen order. The finalize_evidence.py snapshot documents the independent source-count and immutability audit.
