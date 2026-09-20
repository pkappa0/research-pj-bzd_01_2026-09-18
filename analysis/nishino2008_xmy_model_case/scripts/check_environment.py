from pathlib import Path
import sys,json,hashlib,importlib.metadata as md,shutil
ROOT=Path(__file__).resolve().parents[3];R=Path(__file__).resolve().parents[1];S=ROOT/'runs/20260919_0724_10drug_standardized_redocking'
sys.path.insert(0,str(ROOT))
from src import standardized_redocking as sr
old=json.loads((S/'config/software_versions.json').read_text())
diffs={k:[v,md.version(k)] for k,v in old['packages'].items() if md.version(k)!=v}
assert set(diffs)<= {'pip'},diffs
assert sr.sha(sr.VINA)==old['vina_binary_sha256']
base=Path('/opt/homebrew/Caskroom/miniforge/base/envs/plip_env')
current={json.loads(p.read_text())['name']:json.loads(p.read_text()) for p in (base/'conda-meta').glob('*.json')}
condadiff=[]
for a in json.loads((S/'config/conda_base_lock.json').read_text()):
 b=current.get(a['name'],{})
 for k in ['version','build','md5','sha256']:
  if a.get(k) is not None and a.get(k)!=b.get(k):condadiff.append([a['name'],k,a.get(k),b.get(k)])
assert not condadiff,condadiff
scratch=ROOT.parent/'nishino_environment_replay'
for d in ['raw/ligands','raw/receptors','config','tables']:(scratch/d).mkdir(parents=True,exist_ok=True)
shutil.copy2(S/'config/docking_config.json',scratch/'config/docking_config.json')
sr.DRUGS=['diazepam','triazolam'];sr.runpath=lambda:scratch
sr.prepare()
checks=[]
for f in [*(S/'raw/receptors').glob('*'),*(S/'raw/ligands').glob('diazepam.*'),*(S/'raw/ligands').glob('triazolam.*')]:
 dest=scratch/f.relative_to(S);match=sr.sha(f)==sr.sha(dest)
 # SDF header may include timestamp; PDBQT/coordinates are decisive.
 semantic=lambda p:'\n'.join(l for l in p.read_text().splitlines() if not l.startswith(('COMPND','REMARK  Name =')))
 checks.append({'scientific_records_identical':semantic(f)==semantic(dest),'file':str(f.relative_to(S)),'source_sha256':sr.sha(f),'replay_sha256':sr.sha(dest),'exact_match':match})
assert all(c['scientific_records_identical'] for c in checks),checks
for name in ['standardized_redocking.py','standardized_plip.py','standardized_structure_analysis.py']:
 assert (ROOT/'src'/name).read_bytes()==(S/'code'/name).read_bytes(),name
pc=sr.configure_plip();prior=json.loads((S/'config/plip_effective_config.json').read_text())
# Config path values are operational output settings, not interaction thresholds.
params={k:v for k,v in vars(pc).items() if k.isupper() and isinstance(v,(int,float,str,bool,list,dict,type(None)))}
pdiff={k:[v,params.get(k)] for k,v in prior.items() if params.get(k)!=v}
assert not pdiff,pdiff
out={'decision':'reuse frozen diazepam/triazolam; add only brotizolam/lormetazepam','package_differences':diffs,'pip_difference_material':False,'conda_differences':condadiff,'vina_binary_identical':True,'replay_checks':checks,'source_code_identical':True,'plip_effective_config_identical':True,'python':sys.version,'executable':sys.executable}
(R/'qc/environment_check.json').write_text(json.dumps(out,indent=2)+'\n')
for name in ['docking_config.json','software_versions.json','complete_environment_lock.txt','conda_base_lock.json','plip_effective_config.json']:
 shutil.copy2(S/'config'/name,R/'config'/('standardized_'+name))
for name in ['standardized_redocking.py','standardized_plip.py','standardized_structure_analysis.py']:shutil.copy2(ROOT/'src'/name,R/'scripts'/('source_'+name))
print(json.dumps(out,indent=2))
