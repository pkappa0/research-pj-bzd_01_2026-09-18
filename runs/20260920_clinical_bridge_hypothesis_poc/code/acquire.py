from pathlib import Path
import json, hashlib, shutil, subprocess, datetime, time, concurrent.futures, threading
import requests, pandas as pd, numpy as np
ROOT=Path(__file__).resolve().parents[3]; R=Path(__file__).resolve().parents[1]
for d in ['tables','config','raw/responses','raw/lineage','reports','figures']: (R/d).mkdir(parents=True,exist_ok=True)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x): Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False))
D=['diazepam','alprazolam','triazolam','zolpidem','lorazepam','clonazepam','midazolam','temazepam','zopiclone','zaleplon']
S=ROOT/'runs/20260919_0724_10drug_standardized_redocking'; M=ROOT/'runs/20260919_1938_phenotype_fingerprint_poc'
X=S/'tables/fingerprint_10drug_multireceptor.csv'; xf=pd.read_csv(X,index_col=0)
shutil.copy2(X,R/'raw/lineage/structural_fingerprint.csv')
for name in ['phenotype_profile_10drug.csv','phenotype_normalized_within_study.csv','compound_identity_linkage.csv']:
 shutil.copy2(M/'tables'/name,R/'raw/lineage'/name)
shutil.copy2(ROOT/'ANALYSIS_CONTRACT.md',R/'config/ANALYSIS_CONTRACT_source.md')
shutil.copy2('/Users/k-atsumi/.codex/attachments/de15e247-240b-43ea-951f-429ebc0fc71a/pasted-text.txt',R/'config/USER_ANALYSIS_CONTRACT.txt') if Path('/Users/k-atsumi/.codex/attachments/de15e247-240b-43ea-951f-429ebc0fc71a/pasted-text.txt').exists() else None
pd.DataFrame([{'source_run':str(S.relative_to(ROOT)),'source_file':str(X.relative_to(ROOT)),'sha256':sha(X),'feature_count':len(xf.columns),'drug_count':len(xf),'feature_order':json.dumps(xf.columns.tolist()),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}]).to_csv(R/'tables/structural_fingerprint_frozen_manifest.csv',index=False)
axes={'C1':['Ataxia','Balance disorder','Coordination abnormal','Gait disturbance'], 'C2':['Fall'], 'C3':['Somnolence','Sedation'], 'C4':['Dizziness','Vertigo'], 'C5':['Muscular weakness','Hypotonia'], 'C6':['Amnesia','Memory impairment','Cognitive disorder']}
# Exact source term strings, no synonym union. Requested variants kept as explicit unmapped rows.
terms=[(a,t) for a,ts in axes.items() for t in ts]
mapping=[{'phenotype_axis':a,'requested_term':t.lower(),'source_term':t,'mapping_rule':'exact openFDA reactionmeddrapt; availability audited from background; not a licensed MedDRA hierarchy assertion','included':True} for a,t in terms]
for a,t in [('C1','abnormal gait'),('C2','tendency to fall'),('C2','injury secondary to fall'),('C3','drowsiness'),('C5','muscle relaxation')]: mapping.append({'phenotype_axis':a,'requested_term':t,'source_term':None,'mapping_rule':'not independently mapped; no inferred synonym or causal proxy','included':False})
pd.DataFrame(mapping).to_csv(R/'tables/clinical_term_mapping.csv',index=False,na_rep='NA')
config={'axes':axes,'feature_order':[t for _,t in terms],'zero_cells':'if any cell=0 add 0.5 to all four cells; no correction otherwise','minimum_event_reports':5,'unstable_primary_policy':'a<5 remains in raw ROR table but Y primary is NA; all-estimable sensitivity retained','missing':'NA, never impute zero','clinical_similarity':'cosine logROR, pairwise complete fixed terms; require >=7/14 shared terms; report common count','structural_similarity':'cosine all unchanged 48 frequency features; pairwise complete if missing','feature_correlations':'drug-level Spearman, pairwise complete; n<5 insufficient; constant feature rho=NA; no p-value ranking','axis_summary':'none; preserve term-level features','time_window':'receivedate 20040101 through 20251231 inclusive','in_vivo':'archived context-specific values only; no cross-study pooling','pair_statistic':'Spearman over pairs descriptive; 9999 drug-label permutations with fixed seed 20260920 for contextual sensitivity, not n=45 independent p-value','selection':'no clinical-outcome-based feature or term reselection'}
dump(R/'config/clinical_phenotype_config.json',config)
fc={'endpoint':'https://api.fda.gov/drug/event.json','base_query':'receivedate:[20040101 TO 20251231]','drug_query':'patient.drug.openfda.generic_name.exact:UPPERCASE OR patient.drug.medicinalproduct.exact:UPPERCASE','role':'all roles primary fallback; primary suspect not available as PS/SS in openFDA; no invalid cross-drug nested role filtering','deduplication':'API report-document counts; OR union avoids alias double-addition; API revisions/duplicate cases not independently deduplicated','background':'all API report documents in same date range, including reports without harmonized names','demographics':'separate drug-level marginal distributions; not joint drug-event adjusted or report-level microdata','exclusions':'unmapped brands and combinations lacking exact generic identity are not forcibly assigned; eszopiclone is not zopiclone','sensitivity':'harmonized generic-only counts; all-count Y; no suspect-specific ROR claimed'}
dump(R/'config/faers_query_config.json',fc)
dump(R/'STRUCTURAL_RESULTS_FROZEN.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':sha(X),'features':len(xf.columns),'contract_sha256':sha(ROOT/'ANALYSIS_CONTRACT.md'),'config_sha256':sha(R/'config/clinical_phenotype_config.json')})
(R/'config/CONTRACT_ADDENDUM.md').write_text('User-authorized 2026-09-20 extension: primary axis is fixed multi-receptor X versus independently frozen clinical Y. Same-receptor drug contrasts retained; alpha1/alpha2 not replicates. No docking, reselection, ML, mediation or causal inference. Original analysis contract is snapshotted and unchanged. Clinical terms and all analysis rules are fixed before clinical retrieval.\n')
lock=threading.Lock(); provenance=[]
def api(key,params):
 p=R/'raw/responses'/f'{key}.json'; meta=p.with_suffix('.provenance.json')
 if p.exists() and meta.exists():
  j=json.loads(p.read_text()); m=json.loads(meta.read_text())
 else:
  for attempt in range(6):
   try:
    res=requests.get(fc['endpoint'],params=params,timeout=90)
    if res.status_code in [429,500,502,503,504]: time.sleep(2**attempt+1); continue
    j=res.json(); m={'key':key,'url':res.url,'params':params,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'http_status':res.status_code}
    if res.status_code!=200 and not (res.status_code==404 and j.get('error',{}).get('code')=='NOT_FOUND'): raise RuntimeError(str(j))
    dump(p,j); m['sha256']=sha(p); dump(meta,m); break
   except (requests.RequestException,ValueError):
    if attempt==5: raise
    time.sleep(2**attempt)
  else: raise RuntimeError('API failed '+key)
 with lock: provenance.append(m)
 return j,m
def total(key,q):
 j,m=api(key,{'search':q,'limit':1}); return j.get('meta',{}).get('results',{}).get('total',0),m
B=fc['base_query']
def dq(d,generic=False):
 g=f'patient.drug.openfda.generic_name.exact:"{d.upper()}"'
 return g if generic else '('+g+f' OR patient.drug.medicinalproduct.exact:"{d.upper()}")'
def eq(t): return f'patient.reaction.reactionmeddrapt.exact:"{t.upper()}"'
N,nmeta=total('background_total',B)
background={}
for a,t in terms: background[t]=total('background_'+t.replace(' ','_'),B+' AND '+eq(t))[0]
print('Background',N,background,flush=True)
ident=pd.read_csv(S/'tables/ligand_preparation_qc.csv').set_index('drug_id')
# Brand values are evidence-derived from API harmonized field counts, not used as separate query operands.
rows=[]; drugmap=[]; margins=[]
def perdrug(d):
 nd,dm=total(d+'_total',B+' AND '+dq(d)); ng,_=total(d+'_generic_total',B+' AND '+dq(d,True))
 rr=[]; mm=[]; brands=[]
 for field in ['patient.drug.openfda.brand_name.exact','patient.patientsex','patient.patientonsetage','primarysourcecountry.exact','receivedate','serious','patient.drug.drugcharacterization']:
  j,meta=api(d+'_margin_'+field.replace('.','_'),{'search':B+' AND '+dq(d),'count':field,'limit':100})
  for v in j.get('results',[]): mm.append({'drug':d,'field':field,'value':v.get('term',v.get('time')),'report_count':v['count'],'query':meta['url'],'retrieval_date':meta['retrieved_utc'],'scope':'top-100 drug-report marginal (may truncate); drug-array fields include co-reported drugs; age units mixed, not age-years'})
  # Brand count is NOT assigned to target identity because array values include co-medications.
 for a,t in terms:
  n,meta=total(d+'_'+t.replace(' ','_'),B+' AND '+dq(d)+' AND '+eq(t))
  rr.append({'drug':d,'adverse_event_term':t,'phenotype_axis':a,'report_count':n,'total_drug_reports':nd,'background_event_reports':background[t],'background_total':N,'generic_only_total':ng,'role_cod':'NA','suspect_status':'all roles; target-specific suspect not inferred','sex':'ALL; marginal table only','age':'ALL; marginal table only','country':'ALL; marginal table only','report_date_quarter':'2004-01-01 through 2025-12-31 receivedate','seriousness':'ALL; marginal table only','source_dataset':'openFDA FAERS report documents','query':meta['url'],'retrieval_date':meta['retrieved_utc']})
 print(d,nd,'completed',flush=True)
 return rr,mm
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 for rr,mm in ex.map(perdrug,D): rows+=rr; margins+=mm
pd.DataFrame(rows).to_csv(R/'tables/clinical_raw_counts.csv',index=False)
pd.DataFrame(margins).to_csv(R/'tables/clinical_demographic_marginals.csv',index=False)
# ChEMBL synonyms come from original frozen acquisition, if available; otherwise fetch identity endpoint.
for d in D:
 cid=ident.loc[d,'molecule_chembl_id']; url=f'https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json'; p=R/'raw'/f'chembl_{cid}.json'
 try:
  if not p.exists():
   res=requests.get(url,timeout=90); res.raise_for_status(); p.write_text(res.text)
  cj=json.loads(p.read_text()); synonyms=cj.get('molecule_synonyms',[])
  brand=[x['molecule_synonym'] for x in synonyms if x.get('syn_type')=='TRADE_NAME']
  variants=[x['molecule_synonym'] for x in synonyms if x.get('syn_type')!='TRADE_NAME']
 except Exception as e: brand=[]; variants=[]; print('synonym unavailable',d,str(e),flush=True)
 drugmap.append({'drug':d,'generic_name':d,'chembl_id':cid,'standard_inchikey':ident.loc[d,'standard_inchikey'],'brand_synonyms':json.dumps(brand),'spelling_variants_and_other_synonyms':json.dumps(variants),'query_key':dq(d),'synonym_use':'documented only, not independently summed; exact generic or medicinalproduct primary','source':url,'identity_source':str(S.relative_to(ROOT))+'/tables/ligand_preparation_qc.csv'})
pd.DataFrame(drugmap).to_csv(R/'tables/drug_name_mapping.csv',index=False)
dump(R/'raw/query_provenance.json',provenance)
raw=pd.DataFrame(rows); out=[]
for r in rows:
 a=r['report_count']; b=r['total_drug_reports']-a; c=r['background_event_reports']-a; dd=N-a-b-c
 assert min(a,b,c,dd)>=0
 vals=np.array([a,b,c,dd],float); corrected=(vals==0).any()
 if corrected: vals+=.5
 logror=np.log(vals[0])+np.log(vals[3])-np.log(vals[1])-np.log(vals[2]); se=np.sqrt((1/vals).sum())
 usable=r['total_drug_reports']>0 and r['background_event_reports']>0
 out.append({**r,'a':a,'b':b,'c':c,'d':dd,'ROR':np.exp(logror) if usable else np.nan,'logROR':logror if usable else np.nan,'CI95_low':np.exp(logror-1.96*se) if usable else np.nan,'CI95_high':np.exp(logror+1.96*se) if usable else np.nan,'zero_cell_corrected':corrected,'unstable_signal':a<5,'estimable':usable})
out=pd.DataFrame(out);out.to_csv(R/'tables/clinical_disproportionality.csv',index=False,na_rep='NA')
y=out.pivot(index='drug',columns='adverse_event_term',values='logROR').reindex(index=D,columns=config['feature_order'])
y.to_csv(R/'tables/clinical_fingerprint_all_estimable.csv',na_rep='NA')
cnt=out.pivot(index='drug',columns='adverse_event_term',values='a').reindex_like(y);y=y.mask(cnt<5)
y.to_csv(R/'tables/clinical_fingerprint_10drug.csv',na_rep='NA')
qc=pd.DataFrame({'drug':D,'nonmissing_features':y.notna().sum(axis=1).values,'missing_features':y.isna().sum(axis=1).values,'total_drug_reports':[raw[raw.drug==d].total_drug_reports.iloc[0] for d in D],'unstable_cells':[(cnt.loc[d]<5).sum() for d in D]})
qc.to_csv(R/'tables/clinical_missingness_summary.csv',index=False)
pd.DataFrame([{'axis':a,'term':t,'background_reports':background[t],'stable_drugs':int(y[t].notna().sum()),'source_term_available':background[t]>0} for a,t in terms]).to_csv(R/'tables/clinical_term_coverage.csv',index=False)
dump(R/'CLINICAL_RESULTS_FROZEN.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'clinical_sha256':sha(R/'tables/clinical_fingerprint_10drug.csv'),'all_estimable_sha256':sha(R/'tables/clinical_fingerprint_all_estimable.csv'),'config_sha256':sha(R/'config/clinical_phenotype_config.json'),'query_config_sha256':sha(R/'config/faers_query_config.json'),'qc':qc.to_dict('records'),'background_total':N,'term_coverage':background,'note':'Y-only QC complete before correspondence. No terms reselected.'})
print('Y FROZEN',flush=True)
