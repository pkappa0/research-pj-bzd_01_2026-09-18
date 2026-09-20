"""Fixed structural-only extension; never reads M or Y."""
from pathlib import Path
import sys,json,shutil,csv,tarfile
ROOT=Path(__file__).resolve().parents[3]; R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src import standardized_redocking as sr
from src import standardized_plip as sp
from concurrent.futures import ThreadPoolExecutor,as_completed
S=ROOT/'runs/20260919_0724_10drug_standardized_redocking';W=ROOT.parent/'nishino_jobs'
NEW=['brotizolam','lormetazepam'];ALL=['diazepam','triazolam',*NEW]
for mod in [sr,sp]:mod.ROOT=ROOT;mod.WORK=W;mod.DRUGS=NEW;mod.runpath=lambda:R
W.mkdir(exist_ok=True)
def prepare():
 assert json.loads((R/'qc/environment_check.json').read_text())['vina_binary_identical']
 cfg=json.loads((S/'config/docking_config.json').read_text());cfg['drugs']=ALL;cfg['ligand']['source']=str((R/'data/preparation_identity_source.csv').relative_to(ROOT));sr.save(R/'config/docking_config.json',cfg)
 rows=[r for r in sr.read(ROOT/'runs/20260918_1117_chembl_acquisition/processed/compounds.csv') if r['drug_id'] in ALL[:2]]
 rows=[{k:r[k] for k in ['drug_id','molecule_chembl_id','canonical_smiles','standard_inchikey','source_file','source_sha256']} for r in rows]
 for d in NEW:
  p=R/'raw'/f'chembl_{d}.json';ms=json.loads(p.read_text())['molecules'];assert len(ms)==1
  j=ms[0];st=j['molecule_structures'];rows.append(dict(drug_id=d,molecule_chembl_id=j['molecule_chembl_id'],canonical_smiles=st['canonical_smiles'],standard_inchikey=st['standard_inchi_key'],source_file=str(p.relative_to(ROOT)),source_sha256=sr.sha(p)))
 sr.table(R/'data/preparation_identity_source.csv',rows)
 # Exactly the original ligand preparation; receptor replay is replaced with
 # byte-identical frozen prepared receptors after verification (header paths only).
 sr.prepare()
 shutil.copy2(R/'tables/receptor_preparation_qc.csv',R/'tables/receptor_preparation_replay_qc.csv')
 for p in (S/'raw/receptors').iterdir():shutil.copy2(p,R/'raw/receptors'/p.name)
 sr.table(R/'tables/receptor_preparation_qc.csv',[dict(row,production_policy='byte-identical frozen prepared receptor reused') for row in sr.read(S/'tables/receptor_preparation_qc.csv')])
 for d in ALL[:2]:
  for ext in ['sdf','pdbqt']:shutil.copy2(S/'raw/ligands'/f'{d}.{ext}',R/'raw/ligands'/f'{d}.{ext}')
 sr.save(R/'config/PREPARED_INPUTS_FROZEN.json',{'files':{str(p.relative_to(R)):sr.sha(p) for folder in ['raw/receptors','raw/ligands'] for p in (R/folder).iterdir() if p.is_file()},'source_receptors_byte_identical':True})
 manifest=[]
 for row in rows:
  manifest.append(dict(drug=row['drug_id'],canonical_smiles=row['canonical_smiles'],source_database='ChEMBL',source_identifier=row['molecule_chembl_id'],standard_inchikey=row['standard_inchikey'],structure_preparation_method=cfg['ligand']['embedding']+'; '+cfg['ligand']['minimization'],charge_method='Meeko 0.7.1 Gasteiger; as-recorded formal charge; no pH titration',notes=cfg['ligand']['stereochemistry'],source_file=row['source_file'],source_sha256=row['source_sha256']))
 sr.table(R/'data/nishino_primary4_ligand_manifest.csv',manifest)
def dock():sr.docking()
def plip():
 tasks=[(d,b,s) for d in NEW for b in sr.BLOCKS for s in sr.SEEDS]
 def call(t):
  import subprocess
  log=W/'jobs'/t[0]/t[1]/str(t[2])/'plip.log'
  with log.open('w') as f:r=subprocess.run([sys.executable,str(Path(__file__)),'worker',*map(str,t)],stdout=f,stderr=subprocess.STDOUT)
  assert r.returncode==0,str(t)
  return t
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(call,t) for t in tasks]):print('PLIP complete',f.result(),flush=True)
 for drug in NEW:
  with tarfile.open(R/'raw/archives'/f'{drug}.tar.gz','w:gz') as tf:tf.add(W/'jobs'/drug,arcname=drug)
 poses=[];contacts=[]
 for t in tasks:
  j=W/'jobs'/t[0]/t[1]/str(t[2]);poses+=json.loads((j/'pose_manifest.json').read_text());contacts+=json.loads((j/'interactions.json').read_text())
 assert all(r['status']=='success' for r in poses)
 sr.table(R/'tables/pose_manifest_new.csv',poses);sr.table(R/'tables/plip_interactions_new.csv',contacts)
 sr.save(R/'raw/new_contacts.json',contacts);sr.save(R/'raw/new_poses.json',poses)
if __name__=='__main__':
 if sys.argv[1]=='worker':sp.worker(sys.argv[2],sys.argv[3],int(sys.argv[4]))
 else:{'prepare':prepare,'dock':dock,'plip':plip}[sys.argv[1]]()
