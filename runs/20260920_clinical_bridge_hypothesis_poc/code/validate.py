from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 f=json.loads((R/'CLINICAL_RESULTS_FROZEN.json').read_text());s=json.loads((R/'STRUCTURAL_RESULTS_FROZEN.json').read_text())
 assert sha(R/'tables/clinical_fingerprint_10drug.csv')==f['clinical_sha256']
 assert sha(R/'config/clinical_phenotype_config.json')==f['config_sha256']
 assert sha(R/'config/faers_query_config.json')==f['query_config_sha256']
 assert sha(R/'raw/lineage/structural_fingerprint.csv')==s['sha256']
 a=pd.read_csv(R/'tables/clinical_disproportionality.csv');assert len(a)==140
 assert not a.duplicated(['drug','adverse_event_term']).any()
 assert np.all(a.a+a.b==a.total_drug_reports)
 assert np.all(a.a+a.c==a.background_event_reports)
 assert np.all(a.a+a.b+a.c+a.d==a.background_total)
 v=a[['a','b','c','d']].to_numpy(float);v[(v==0).any(axis=1)]+=.5
 assert np.allclose(a.ROR,(v[:,0]*v[:,3])/(v[:,1]*v[:,2]))
 assert np.all(a.CI95_low<=a.ROR) and np.all(a.CI95_high>=a.ROR)
 y=pd.read_csv(R/'tables/clinical_fingerprint_10drug.csv',index_col=0)
 for row in a.itertuples():assert pd.isna(y.loc[row.drug,row.adverse_event_term]) == (row.a<5 or not row.estimable)
 p=pd.read_csv(R/'tables/structural_clinical_pairwise_similarity.csv');assert len(p)==45
 assert len({tuple(sorted([r.drug_1,r.drug_2])) for r in p.itertuples()})==45
 c=pd.read_csv(R/'tables/structural_feature_clinical_correlations.csv');assert len(c)==48*14
 assert c[c.n_effective<5].rho.isna().all()
 assert c.structural_feature.nunique()==48
 assert f['created_utc']>=s['created_utc']
 for v in json.loads((R/'raw/query_provenance.json').read_text()):
  assert sha(R/'raw/responses'/f"{v['key']}.json")==v['sha256']
 print('PASS: frozen hashes, prespecified config, 140 independent contingency reconstructions, CIs, NA policy, 45 unique pairs, all 672 correlations and raw response integrity')
if __name__=='__main__':main()
