"""Fetch ChEMBL molecule and GABA-A-related activity data.

The input drug list is authoritative. Fuzzy molecule matches are retained for QC
but never selected automatically. Full API responses are archived as JSONL/pages.
"""
from __future__ import annotations
import argparse, csv, json, re, time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, urljoin, urlparse
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import pandas as pd
import yaml

MOLECULE_COLUMNS = ['drug_id','input_drug_name','molecule_chembl_id','preferred_name','canonical_smiles','standard_inchikey','molecule_type','chirality','parent_chembl_id','active_chembl_id','salt_or_multicomponent','stereo_note','match_status','match_method','candidate_count','api_retrieved_at','chembl_version']
ACTIVITY_COLUMNS = ['drug_id','input_drug_name','molecule_chembl_id','preferred_name','canonical_smiles','standard_inchikey','target_chembl_id','target_name','target_type','target_class','receptor_subtype','assay_chembl_id','assay_type','standard_type','standard_relation','standard_value','standard_units','pchembl_value','organism','document_chembl_id','reference','source','src_id','activity_id','record_id','api_retrieved_at','chembl_version']
QC_COLUMNS = ['drug_id','input_drug_name','status','match_status','match_method','molecule_candidate_count','molecule_chembl_id','activity_count','target_candidate_count','notes']

def now(): return datetime.now(timezone.utc).isoformat()
def norm(s): return re.sub(r'\s+', ' ', str(s or '').strip()).casefold()

class ChemblClient:
    def __init__(self, base, timeout=45, delay=0.1, raw_dir=None):
        self.base=base.rstrip('/')+'/'
        self.timeout=float(timeout); self.delay=float(delay); self.raw_dir=Path(raw_dir) if raw_dir else None; self.page_no={}
    def get(self, endpoint, params=None, raw_label=None):
        if endpoint.startswith('http'): url=endpoint
        elif endpoint.startswith('/chembl/api/'):
            parsed=urlparse(self.base); url=f'{parsed.scheme}://{parsed.netloc}{endpoint}'
        else: url=urljoin(self.base, endpoint.lstrip('/'))
        if params: url += ('&' if '?' in url else '?') + urlencode(params)
        req=Request(url, headers={'Accept':'application/json','User-Agent':'gaba-a-fingerprint-poc/0.1'})
        last_exc=None
        for attempt in range(3):
            try:
                with urlopen(req, timeout=self.timeout) as r: data=json.load(r)
                break
            except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
                last_exc=exc
                if attempt < 2: time.sleep(2 ** attempt)
        else:
            raise RuntimeError(f'ChEMBL request failed: {url}: {last_exc}') from last_exc
        if raw_label and self.raw_dir:
            d=self.raw_dir/'api_pages'; d.mkdir(parents=True, exist_ok=True); i=self.page_no.get(raw_label,0); self.page_no[raw_label]=i+1
            (d/f'{raw_label}_{i:04d}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2))
        if self.delay: time.sleep(self.delay)
        return data
    def all(self, endpoint, params=None, raw_label=None):
        params=dict(params or {}); params.setdefault('limit',1000); params.setdefault('offset',0); records=[]; first=True
        while True:
            data=self.get(endpoint, params if first else None, raw_label=raw_label); first=False
            key=next((k for k in data if k not in ('page_meta','_meta')), None)
            if key and isinstance(data[key],list): records.extend(data[key])
            nxt=(data.get('page_meta') or {}).get('next')
            if not nxt: break
            endpoint=nxt; params=None
        return records

def target_is_gabaa(t):
    name=norm(t.get('pref_name'))
    if 'associated protein' in name or 'rho' in name or 'gaba-c' in name: return False
    return ('gaba-a' in name or 'gaba a receptor' in name or 'gaba(a)' in name or 'gamma-aminobutyric acid receptor' in name)
def subtype_from_name(name):
    text=norm(name)
    if 'subunit' in text: text=text.split('subunit',1)[1]
    elif ';' in text: text=text.split(';',1)[1]
    vals=[]
    for family,num in re.findall(r'\b(alpha|beta|gamma|delta|epsilon|theta|pi)(?:[- ]?(\d+))?', text):
        vals.append(f'{family}{num}')
    return '/'.join(dict.fromkeys(vals))
def classify_target(t):
    name=norm(t.get('pref_name')); typ=str(t.get('target_type') or '')
    if not subtype_from_name(name) and ('gaba-a receptor' in name or 'gaba a receptor' in name or 'gaba(a)' in name): return 'generic_gabaa_target'
    if '/' in name or 'complex' in name or 'interface' in name or 'channel' in name or 'complex' in typ.casefold(): return 'composite_target'
    if 'subunit' in name: return 'subunit_specific'
    return 'other_gabaa_related'
def molecule_note(m):
    smi=((m.get('molecule_structures') or {}).get('canonical_smiles') or '')
    notes=[]
    if '.' in smi: notes.append('multi_component_or_salt_like_structure')
    if m.get('chirality') not in (None,0,'0'): notes.append(f"chirality={m.get('chirality')}")
    h=m.get('molecule_hierarchy') or {}
    if h.get('parent_chembl_id') and h.get('parent_chembl_id') != m.get('molecule_chembl_id'): notes.append(f"parent={h.get('parent_chembl_id')}")
    return ';'.join(notes)

def read_drug_list(path):
    p=Path(path)
    if not p.exists(): return None, 'input drug list is missing'
    d=pd.read_csv(p, dtype='string')
    missing=[c for c in ['drug_id','input_drug_name'] if c not in d.columns]
    if missing: return None, 'missing required columns: '+','.join(missing)
    d=d[['drug_id','input_drug_name']].copy(); d['drug_id']=d.drug_id.str.strip().str.lower(); d['input_drug_name']=d.input_drug_name.str.strip(); return d, ''

def resolve_molecule(client, name):
    exact=client.all('/molecule.json', {'pref_name__iexact':name}, raw_label='molecule_exact')
    if len(exact)==1: return exact[0], exact, 'exact_preferred_name'
    candidates=client.all('/molecule/search.json', {'q':name}, raw_label='molecule_search')
    n=norm(name); exact_syn=[]
    for m in candidates:
        if norm(m.get('pref_name'))==n: exact_syn.append(m)
        for syn in m.get('molecule_synonyms') or []:
            if norm(syn.get('synonyms') or syn.get('molecule_synonym'))==n: exact_syn.append(m); break
    uniq={m.get('molecule_chembl_id'):m for m in exact_syn if m.get('molecule_chembl_id')}
    if len(uniq)==1: return next(iter(uniq.values())), candidates, 'exact_synonym'
    return None, candidates, 'ambiguous_exact' if len(uniq)>1 else 'unresolved_fuzzy'

def fetch_activity_batches(client, molecule_id, target_ids, drug_id):
    """Fetch in bounded target batches; very large activity pages can time out."""
    records=[]; errors=[]
    for batch_no in range(0, len(target_ids), 20):
        batch=target_ids[batch_no:batch_no+20]
        try:
            records.extend(client.all('/activity.json', {
                'molecule_chembl_id': molecule_id,
                'target_chembl_id__in': ','.join(batch),
                'limit': 10,
            }, raw_label=f'activity_{drug_id}_chunk{batch_no//20:02d}'))
        except RuntimeError as exc:
            errors.append(f'batch_{batch_no//20:02d}: {exc}')
    return records, errors

def fetch(cfg):
    ccfg=cfg['chembl']; raw=Path(ccfg['raw_output_dir']); raw.mkdir(parents=True,exist_ok=True); retrieved=now(); client=ChemblClient(ccfg['api_base'],ccfg.get('request_timeout_seconds',45),ccfg.get('request_delay_seconds',0.1),raw)
    metadata={'api_base':ccfg['api_base'],'retrieved_at_utc':retrieved,'chembl_version':None,'status':'unknown','target_search_terms':ccfg.get('target_search_terms',[])}
    try:
        status=client.get('/status.json',raw_label='status'); metadata['chembl_version']=status.get('chembl_db_version'); metadata['chembl_release_date']=status.get('chembl_release_date'); metadata['status']=status.get('status')
    except RuntimeError as exc: metadata['status']='api_error'; metadata['error']=str(exc)
    (raw/'chembl_api_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2))
    targets=[]; target_errors=[]
    for term in ccfg.get('target_search_terms',[]):
        try: targets.extend(client.all('/target.json',{'pref_name__icontains':term},raw_label='target_search'))
        except RuntimeError as exc: target_errors.append(str(exc))
    for tid in ccfg.get('target_chembl_ids',[]):
        try: targets.append(client.get(f'/target/{tid}.json',raw_label='target_detail'))
        except RuntimeError as exc: target_errors.append(str(exc))
    target_map={t.get('target_chembl_id'):t for t in targets if t.get('target_chembl_id') and target_is_gabaa(t)}
    Path(ccfg['target_candidates_output']).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{'target_chembl_id':t.get('target_chembl_id'),'target_name':t.get('pref_name'),'target_type':t.get('target_type'),'organism':t.get('organism'),'receptor_subtype':subtype_from_name(t.get('pref_name')),'target_class':classify_target(t)} for t in target_map.values()]).to_csv(ccfg['target_candidates_output'],index=False)
    drug_list, list_error=read_drug_list(ccfg['input_drug_list']); molecules=[]; activities=[]; qc=[]; source_cache={}; raw_jsonl=(raw/'molecule_raw.jsonl').open('w'); act_jsonl=(raw/'activity_raw.jsonl').open('w')
    if drug_list is None:
        empty=pd.DataFrame(columns=ACTIVITY_COLUMNS)
        pd.DataFrame(columns=MOLECULE_COLUMNS).to_csv(raw/'molecules.csv',index=False); empty.to_csv(raw/'activities.csv',index=False); empty.to_csv(ccfg['processed_output'],index=False); pd.DataFrame([{'drug_id':'','input_drug_name':'','status':'missing_input','match_status':'missing_input','match_method':'','molecule_candidate_count':0,'molecule_chembl_id':'','activity_count':0,'target_candidate_count':len(target_map),'notes':list_error}]).to_csv(ccfg['qc_output'],index=False); raw_jsonl.close(); act_jsonl.close(); return {'molecules':0,'activities':0,'targets':len(target_map),'status':'missing_input'}
    for _, row in drug_list.iterrows():
        did=str(row.drug_id); name=str(row.input_drug_name); selected,candidates,method=resolve_molecule(client,name); candidate_ids='|'.join(str(m.get('molecule_chembl_id')) for m in candidates if m.get('molecule_chembl_id'))
        for c in candidates: raw_jsonl.write(json.dumps({'input_drug_name':name,'drug_id':did,'search_candidate':c},ensure_ascii=False)+'\n')
        if selected is None:
            qc.append({'drug_id':did,'input_drug_name':name,'status':'unresolved_name','match_status':method,'match_method':method,'molecule_candidate_count':len(candidates),'molecule_chembl_id':'','activity_count':0,'target_candidate_count':len(target_map),'notes':'No unique exact preferred-name/synonym match; candidates='+candidate_ids}); continue
        mid=selected['molecule_chembl_id']; detail=client.get(f'/molecule/{mid}.json',raw_label='molecule_detail'); sm=detail.get('molecule_structures') or {}; hier=detail.get('molecule_hierarchy') or {}; raw_jsonl.write(json.dumps({'input_drug_name':name,'drug_id':did,'selected_molecule':detail},ensure_ascii=False)+'\n'); mrow={'drug_id':did,'input_drug_name':name,'molecule_chembl_id':mid,'preferred_name':detail.get('pref_name'),'canonical_smiles':sm.get('canonical_smiles'),'standard_inchikey':sm.get('standard_inchi_key'),'molecule_type':detail.get('molecule_type'),'chirality':detail.get('chirality'),'parent_chembl_id':hier.get('parent_chembl_id'),'active_chembl_id':hier.get('active_chembl_id'),'salt_or_multicomponent':bool('.' in (sm.get('canonical_smiles') or '')),'stereo_note':molecule_note(detail),'match_status':'matched','match_method':method,'candidate_count':len(candidates),'api_retrieved_at':retrieved,'chembl_version':metadata.get('chembl_version')}; molecules.append(mrow)
        target_ids=sorted(target_map, key=lambda tid: (0 if ('gaba-a receptor;' in norm(target_map[tid].get('pref_name')) or 'gaba a receptor' in norm(target_map[tid].get('pref_name'))) else 1, tid))
        acts, activity_errors=fetch_activity_batches(client,mid,target_ids,did)
        target_errors.extend([f'{did}: {e}' for e in activity_errors])
        for a in acts:
            act_jsonl.write(json.dumps({'drug_id':did,'input_drug_name':name,'activity':a},ensure_ascii=False)+'\n'); t=target_map.get(a.get('target_chembl_id'),{}); ref=' '.join(str(x) for x in [a.get('document_chembl_id'),a.get('document_journal'),a.get('document_year')] if x not in (None,'')); sid=a.get('src_id');
            if sid and sid not in source_cache:
                try: source_cache[sid]=client.get(f'/source/{sid}.json',raw_label='source_detail')
                except RuntimeError: source_cache[sid]={}
            source_obj=source_cache.get(sid,{}) or {}; source_name=source_obj.get('src_short_name') or source_obj.get('src_description') or str(sid or '')
            activities.append({'drug_id':did,'input_drug_name':name,'molecule_chembl_id':mid,'preferred_name':detail.get('pref_name'),'canonical_smiles':sm.get('canonical_smiles'),'standard_inchikey':sm.get('standard_inchi_key'),'target_chembl_id':a.get('target_chembl_id'),'target_name':a.get('target_pref_name') or t.get('pref_name'),'target_type':t.get('target_type'),'target_class':classify_target(t),'receptor_subtype':subtype_from_name(a.get('target_pref_name') or t.get('pref_name')),'assay_chembl_id':a.get('assay_chembl_id'),'assay_type':a.get('assay_type'),'standard_type':a.get('standard_type'),'standard_relation':a.get('standard_relation'),'standard_value':a.get('standard_value'),'standard_units':a.get('standard_units'),'pchembl_value':a.get('pchembl_value'),'organism':a.get('target_organism'),'document_chembl_id':a.get('document_chembl_id'),'reference':ref,'source':source_name,'src_id':sid,'activity_id':a.get('activity_id'),'record_id':a.get('record_id'),'api_retrieved_at':retrieved,'chembl_version':metadata.get('chembl_version')})
        note='No GABA-A-related activity found' if not acts else ''
        if activity_errors: note=(note+'; ' if note else '')+'Some target batches failed; see fetch_errors.json'
        qc.append({'drug_id':did,'input_drug_name':name,'status':'loaded','match_status':'matched','match_method':method,'molecule_candidate_count':len(candidates),'molecule_chembl_id':mid,'activity_count':len(acts),'target_candidate_count':len(target_map),'notes':note})
    raw_jsonl.close(); act_jsonl.close(); pd.DataFrame(molecules,columns=MOLECULE_COLUMNS).to_csv(raw/'molecules.csv',index=False); pd.DataFrame(activities,columns=ACTIVITY_COLUMNS).to_csv(raw/'activities.csv',index=False); processed=pd.DataFrame(activities,columns=ACTIVITY_COLUMNS); processed.to_csv(ccfg['processed_output'],index=False); pd.DataFrame(qc,columns=QC_COLUMNS).to_csv(ccfg['qc_output'],index=False); (raw/'fetch_errors.json').write_text(json.dumps(target_errors,ensure_ascii=False,indent=2))
    return {'molecules':len(molecules),'activities':len(activities),'targets':len(target_map)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',default='config/poc.yaml'); args=ap.parse_args(); cfg=yaml.safe_load(Path(args.config).read_text()); print(json.dumps(fetch(cfg),ensure_ascii=False,indent=2,default=str))

if __name__ == '__main__':
    main()
