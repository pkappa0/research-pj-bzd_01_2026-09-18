# PHARMACOLOGY_INTEGRATION_REPORT

## Scope and rules

ChEMBL raw activity is retained as-is. GtoPdb rows are separate evidence records, and duplicate candidates are flagged rather than silently removed. Ki, Kd, IC50 and EC50 remain distinct; units, species and beta subtypes are not merged.

- Drugs with alpha1/alpha2 quantitative data in GtoPdb alone: **0**
- Drugs with alpha1/alpha2 quantitative data in the combined master: **5**
- Tier 1: **4**
- Tier 2 or better: **4**
- Tier 3 or better: **5**
- ChEMBL-strict missing but GtoPdb-complemented: **none in current run**
- GtoPdb without subtype-specific quantitative data: **diazepam, lorazepam, clonazepam, alprazolam, midazolam, temazepam, triazolam, zolpidem, zopiclone, zaleplon**
- Duplicate candidate rows: **0**

The tier is a data comparability classification, not a potency or efficacy ranking.

## Per-drug QC

   drug_id input_drug_name  gtopdb_alpha1_quantitative  gtopdb_alpha2_quantitative  gtopdb_both  combined_alpha1_quantitative  combined_alpha2_quantitative  combined_alpha1_alpha2  best_tier  tier1_pairs  tier2_pairs  tier3_pairs  chEMBL_strict_comparable  gtopdb_subtype_specific_missing
  diazepam        Diazepam                       False                       False        False                          True                          True                    True        1.0           13           83           82                      True                             True
 lorazepam       Lorazepam                       False                       False        False                          True                         False                   False        NaN            0            0            0                     False                             True
clonazepam      Clonazepam                       False                       False        False                          True                         False                   False        NaN            0            0            0                     False                             True
alprazolam      Alprazolam                       False                       False        False                          True                          True                    True        1.0            1            0            2                     False                             True
 midazolam       Midazolam                       False                       False        False                          True                         False                   False        NaN            0            0            0                     False                             True
 temazepam       Temazepam                       False                       False        False                          True                         False                   False        NaN            0            0            0                     False                             True
 triazolam       Triazolam                       False                       False        False                          True                          True                    True        1.0            1            0            0                      True                             True
  zolpidem        Zolpidem                       False                       False        False                          True                          True                    True        1.0            9           45           40                      True                             True
 zopiclone       Zopiclone                       False                       False        False                          True                         False                   False        NaN            0            0            0                     False                             True
  zaleplon        Zaleplon                       False                       False        False                          True                          True                    True        3.0            0            0            1                     False                             True

## Current limitation

If the API key was unavailable and the fallback CSV was absent, these counts describe the existing ChEMBL master only; they do not claim that GtoPdb has no corresponding data. Rerun after providing the key or official CSV.
