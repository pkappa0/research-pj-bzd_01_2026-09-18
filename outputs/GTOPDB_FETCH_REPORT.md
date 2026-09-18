# GTOPDB_FETCH_REPORT

This is an independent GtoPdb supplementation phase. Existing ChEMBL and strict 3-drug outputs were not modified.

- API key status: **absent**
- Fallback CSV status: **missing**
- GtoPdb normalized GABA-A rows: **0**
- GtoPdb version: **not exposed by response**

The GtoPdb REST API was not called because `GTP_API_KEY` is absent or authentication failed. Set the key only in the process environment:

```bash
export GTP_API_KEY='your-key'
python fetch_gtopdb.py --config config/poc.yaml
```

Do not place the key in code, YAML, CSV, reports, or Git history. Alternatively place the official download at `data/raw/gtopdb/gtopdb_interactions.csv` and rerun.
