"""Apply the existing exact-name, all-role, fixed-window clinical method."""
from pathlib import Path
import requests,json,time,datetime,hashlib
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];ROOT=R.parents[1];OLD=ROOT/'runs/20260920_clinical_bridge_hypothesis_poc'
D=['diazepam','triazolam','brotizolam','lormetazepam'];T=['Ataxia','Balance disorder','Coordination abnormal','Gait disturbance']
B='receivedate:[20040101 TO 20251231]';P=R/'raw/openfda';P.mkdir(exist_ok=True)
for name in ['clinical_phenotype_config.json','faers_query_config.json']:(R/'config'/name).write_bytes((OLD/'config'/name).read_bytes())
def get(key,q):
 p=P/(key+'.json');meta=P/(key+'.provenance.json')
 if not p.exists():
  for i in range(6):
   res=requests.get('https://api.fda.gov/drug/event.json',params={'search':q,'limit':1},timeout=60)
   if res.status_code in [429,500,502,503,504]:time.sleep(2**i);continue
   j=res.json();assert res.status_code==200 or (res.status_code==404 and j.get('error',{}).get('code')=='NOT_FOUND'),(res.status_code,j)
   # Retain counts + API metadata only; do not archive incidental patient reports.
   out={'meta':j.get('meta',{}),'error':j.get('error'),'report_total':j.get('meta',{}).get('results',{}).get('total',0)}
   p.write_text(json.dumps(out,indent=2));meta.write_text(json.dumps({'url':res.url,'query':q,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'http_status':res.status_code,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'results_omitted':'counts-only extraction; individual report not retained'},indent=2));break
  else:raise RuntimeError(key)
 return json.loads(p.read_text())['report_total']
def dq(d):return '(patient.drug.openfda.generic_name.exact:"'+d.upper()+'" OR patient.drug.medicinalproduct.exact:"'+d.upper()+'")'
def eq(t):return 'patient.reaction.reactionmeddrapt.exact:"'+t.upper()+'"'
N=get('background_total',B);backs={t:get('background_'+t.replace(' ','_'),B+' AND '+eq(t)) for t in T}
frozen=json.loads((OLD/'CLINICAL_RESULTS_FROZEN.json').read_text());same=(N==frozen['background_total'] and all(backs[t]==frozen['term_coverage'][t] for t in T))
assert same,'API background changed; do not mix snapshots silently'
orig=pd.read_csv(OLD/'tables/clinical_disproportionality.csv',float_precision='round_trip');oldY=pd.read_csv(OLD/'tables/clinical_fingerprint_10drug.csv',float_precision='round_trip').set_index('drug');rows=[]
for d in D:
 nd=get(d+'_total',B+' AND '+dq(d))
 for t in T:
  a=get(d+'_'+t.replace(' ','_'),B+' AND '+dq(d)+' AND '+eq(t));b=nd-a;c=backs[t]-a;dd=N-a-b-c;assert min(a,b,c,dd)>=0
  if d in D[:2]:
   o=orig[(orig.drug==d)&(orig.adverse_event_term==t)].iloc[0];assert (a,b,c,dd)==tuple(int(o[k]) for k in ['a','b','c','d']),'Frozen count mismatch'
   raw=float(o.logROR);value=oldY.loc[d,t];reused=True
  else:
   cells=np.array([a,b,c,dd],float)
   if (cells==0).any():cells+=.5
   raw=float(np.log(cells[0])+np.log(cells[3])-np.log(cells[1])-np.log(cells[2])) if nd>0 and backs[t]>0 else np.nan
   value=raw if a>=5 else np.nan;reused=False
  rows.append(dict(drug=d,clinical_term=t,exposed_target_reports=nd,a_target_event=a,b_target_nonevent=b,c_comparator_event=c,d_comparator_nonevent=dd,comparator_reports=N-nd,background_event_reports=backs[t],total_reports=N,raw_estimable=bool(nd>0 and backs[t]>0),estimable=bool(np.isfinite(value)),missingness_reason='' if np.isfinite(value) else ('no_target_reports' if nd==0 else 'target_event_count_below_5'),logROR=value,raw_logROR_before_low_count_mask=raw,exact_drug_normalization_rule=dq(d),synonym_rule='uppercase exact generic OR medicinalproduct; no brand expansion; no summed synonyms',frozen_Y_reused=reused,date_window='20040101-20251231',roles='all',source='openFDA FAERS report-document counts'))
 print(d,nd,[(r['clinical_term'],r['a_target_event'],r['estimable']) for r in rows if r['drug']==d],flush=True)
out=pd.DataFrame(rows);out.to_csv(R/'qc/nishino_primary4_clinical_coverage.csv',index=False,na_rep='NA')
y=out.pivot(index='drug',columns='clinical_term',values='logROR').loc[D,T];y.columns=['logROR_ataxia','logROR_balance_disorder','logROR_coordination_abnormal','logROR_gait_disturbance'];y.to_csv(R/'data/nishino_primary4_clinical_Y.csv',na_rep='NA')
(R/'qc/clinical_snapshot_check.json').write_text(json.dumps({'background_matches_frozen':same,'existing_two_drug_counts_match':True,'background_total':N,'nonmissing_C1_cells_by_drug':y.notna().sum(axis=1).to_dict(),'no_database_switch':True},indent=2)+'\n')
