"""Frozen 48-feature schema applied without inspecting M or Y."""
from pathlib import Path
import sys,json,shutil,tarfile,importlib.util,datetime,hashlib
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];ROOT=R.parents[1];sys.path.insert(0,str(ROOT))
from src import standardized_structure_analysis as sa
S=ROOT/'runs/20260919_0724_10drug_standardized_redocking';V=ROOT/'runs/20260920_vina_score_layer'
D=['diazepam','triazolam','brotizolam','lormetazepam'];B=['alpha1','alpha2'];SEEDS=[2026091901+i for i in range(5)]
key=pd.read_csv(ROOT/'runs/20260920_wobbling_focused_correspondence/tables/feature_display_key.csv');features=key.feature_label.tolist();assert len(features)==48
key.to_csv(R/'data/feature_display_key.csv',index=False)
for name in ['diazepam','triazolam']:shutil.copy2(S/'raw/archives'/f'{name}.tar.gz',R/'raw/archives'/f'{name}.tar.gz')
shutil.copy2(V/'data/frozen_plif_fingerprint.csv',R/'data/frozen_10drug_plif_source.csv')
idx=sa.mapping_index(R)
newcontacts=json.loads((R/'raw/new_contacts.json').read_text());poses=json.loads((R/'raw/new_poses.json').read_text())
assert all(p['status']=='success' for p in poses)
old=pd.read_csv(R/'data/frozen_10drug_plif_source.csv').set_index('drug_id');x=old.loc[D[:2],features].copy();unmapped=[];outschema=[];counts=[]
for c in newcontacts:
 k=(c['receptor'],c['residue_chain'],int(c['residue_number']))
 if k not in idx:unmapped.append(c);continue
 cp,res=idx[k];c['frozen_feature_label']=f'{c["receptor"]}|{cp}|{res}|any_contact_frequency'
 if c['frozen_feature_label'] not in features:outschema.append(c)
for drug in D[2:]:
 row=[]
 for f in features:
  b=f.split('|')[0];eligible={(p['seed'],p['pose_rank']) for p in poses if p['drug_id']==drug and p['receptor']==b}
  hits={(c['seed'],c['pose_rank']) for c in newcontacts if c['drug_id']==drug and c.get('frozen_feature_label')==f}&eligible
  assert eligible
  row.append(len(hits)/len(eligible));counts.append(dict(drug=drug,feature_label=f,n_contact_poses=len(hits),n_evaluable_poses=len(eligible),contact_frequency=row[-1]))
 x.loc[drug]=row
x.index.name='drug';x.to_csv(R/'data/nishino_primary4_plif.csv')
pd.DataFrame(counts).to_csv(R/'qc/new_plif_denominators.csv',index=False)
for name,rows in [('unmapped_contacts',unmapped),('mapped_contacts_outside_frozen48',outschema)]:
 pd.DataFrame(rows,columns=list(rows[0]) if rows else ['drug_id','receptor','residue_chain','residue_number','residue_name','interaction_type']).to_csv(R/'qc'/(name+'.csv'),index=False)
parserpath=ROOT/'runs/20260920_archived_vina_score_recovery/scripts/recover_vina_scores.py'
shutil.copy2(parserpath,R/'scripts/source_vina_parser.py')
spec=importlib.util.spec_from_file_location('vp',parserpath);vp=importlib.util.module_from_spec(spec);spec.loader.exec_module(vp)
allrows=[];files=[]
oldmanifest=pd.read_csv(S/'tables/docking_run_manifest.csv');newmanifest=pd.read_csv(R/'tables/docking_run_manifest.csv');dm=pd.concat([oldmanifest[oldmanifest.drug_id.isin(D[:2])],newmanifest])
for drug in D:
 archive=R/'raw/archives'/f'{drug}.tar.gz'
 with tarfile.open(archive) as tf:
  for m in tf.getmembers():
   if not m.name.endswith('/poses.pdbqt'):continue
   d,b,s,_=m.name.split('/');s=int(s);payload=tf.extractfile(m).read();src='raw/archives/'+archive.name+'::'+m.name
   rows,issues,n=vp.parse_pdbqt(payload.decode(),d,b,s,src);assert not issues,(src,issues)
   match=dm[(dm.drug_id==d)&(dm.receptor==b)&(dm.seed==s)].iloc[0];h=hashlib.sha256(payload).hexdigest();assert h==match.output_sha256 and n==match.n_poses
   for row in rows:row['source_pdbqt']=src
   allrows+=rows;files.append(dict(drug=d,receptor=b,seed=s,n_poses=n,sha256=h,source_pdbqt=src,reused=drug in D[:2]))
a=pd.DataFrame(allrows);a.to_csv(R/'data/nishino_primary4_vina_all_poses.csv',index=False)
best=a.sort_values(['vina_score_kcal_mol','pose_rank']).groupby(['drug','receptor_block','seed'],sort=False).head(1).sort_values(['drug','receptor_block','seed']);best.to_csv(R/'data/nishino_primary4_vina_seed_best.csv',index=False)
summary=[]
for d in D:
 for b in B:
  v=best[(best.drug==d)&(best.receptor_block==b)].vina_score_kcal_mol;assert len(v)==5
  summary.append(dict(drug=d,receptor_block=b,n_seeds=len(v),median_best_seed_score=v.median(),mean_best_seed_score=v.mean(),sd_best_seed_score=v.std(ddof=1),IQR_best_seed_score=v.quantile(.75)-v.quantile(.25),min_best_seed_score=v.min(),max_best_seed_score=v.max(),vina_favorability=-v.median()))
s=pd.DataFrame(summary);s.to_csv(R/'data/nishino_primary4_vina_summary.csv',index=False);pd.DataFrame(files).to_csv(R/'qc/docking_file_manifest.csv',index=False)
rows=[]
for d in D:
 r={'drug':d}
 for b in B:
  q=s[(s.drug==d)&(s.receptor_block==b)].iloc[0];r.update({b+'_PLIF_reference':f'data/nishino_primary4_plif.csv#drug={d};columns={b}|*',b+'_median_Vina_score':q.median_best_seed_score,b+'_Vina_favorability':q.vina_favorability})
 rows.append(r)
pd.DataFrame(rows).to_csv(R/'data/nishino_primary4_structural_multiview.csv',index=False)
assert np.array_equal(x.loc[D[:2]].values,old.loc[D[:2],features].values)
freeze={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'M_Y_used_to_construct_X':False,'features':features,'files':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'data').glob('nishino_primary4*') if 'clinical' not in p.name},'reused_drugs':D[:2],'new_drugs':D[2:],'n_docking_files':len(files),'n_poses':len(a),'new_PLIP_successful_poses':len(poses),'unmapped_contact_records':len(unmapped),'mapped_contact_records_outside_frozen48':len(outschema),'out_of_schema_policy':'flag and preserve; never add features; report counts are interactions, not unique residues'}
(R/'config/STRUCTURAL_RESULTS_FROZEN.json').write_text(json.dumps(freeze,indent=2)+'\n');print(json.dumps({k:v for k,v in freeze.items() if k not in ['files','features']},indent=2))
