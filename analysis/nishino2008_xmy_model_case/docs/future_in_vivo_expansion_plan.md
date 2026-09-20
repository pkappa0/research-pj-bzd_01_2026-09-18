# Future in vivo expansion plan

Keep one raw observation per drug × study × assay context × dose × time × endpoint, linked to an immutable source location. Required metadata: canonical parent/active-moiety identity and identifiers; species; strain; sex; age/weight where available; route; formulation/vehicle; dose and unit; observation time; assay protocol; rotarod speed/acceleration and rod geometry; training/selection; positive-event or latency definition; n and repeated-measure design; effect value; uncertainty/CI type and level; reported summary-model method; source page/table and extraction QC.

Explicitly distinguish administered prodrug from active metabolites. Preserve reported units and raw counts; do not pool latency with binary failure or equate different dose/time windows. Record absent metadata as unavailable. Include source anomaly flags and double-check nonmonotonic values without silently correcting them.

- Tier 1: multiple drugs in one study under the same experimental context. Priority acquisition tier. Keep study as a context block and compare drugs within it.
- Tier 2: closely matched protocols across studies with compatible species, route, speed, training, endpoint and time. Require a documented comparability audit and heterogeneity assessment before any combined model.
- Tier 3: heterogeneous protocols or insufficient metadata. Retain as qualitative/contextual evidence; do not pool numerical values as equivalent observations.

Nishino demonstrates why reported ED50 alone is insufficient metadata. Collect the dose/time grid, reported probit inputs/time basis and interval definition alongside it. The primary selection criterion for new studies is context completeness and drug overlap, not favorable or positive correlations. No additional animal experiments are performed or automatically proposed as validated designs here.
