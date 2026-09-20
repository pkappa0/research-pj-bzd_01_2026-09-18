from pathlib import Path
import requests,json,hashlib,datetime
R=Path(__file__).resolve().parents[1]
for name in ['brotizolam','lormetazepam','rilmazafone']:
 url='https://www.ebi.ac.uk/chembl/api/data/molecule.json';p=R/'raw'/f'chembl_{name}.json'
 if not p.exists():
  res=requests.get(url,params={'pref_name__iexact':name,'limit':10},timeout=90);res.raise_for_status();p.write_text(res.text)
  p.with_suffix('.provenance.json').write_text(json.dumps({'url':res.url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
 j=json.loads(p.read_text());print(name,[(x['molecule_chembl_id'],x['pref_name'],x['molecule_structures']['canonical_smiles']) for x in j['molecules']])
