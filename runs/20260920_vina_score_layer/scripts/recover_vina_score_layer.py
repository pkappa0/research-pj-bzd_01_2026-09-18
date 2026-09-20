#!/usr/bin/env python3
"""Read archived PDBQT bytes, never execute docking or recalculate PLIF."""
from pathlib import Path
import hashlib,json,re,shutil,tarfile,math,datetime
from scipy.stats import spearmanr
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
RUN=Path(__file__).resolve().parents[1];ROOT=RUN.parents[1]
SOURCE=ROOT/'runs/20260919_0724_10drug_standardized_redocking'
ORDER_SOURCE=ROOT/'runs/20260920_wobbling_focused_correspondence/tables/drug_display_order.csv'
for d in ['data','data/pdbqt','config','figures']: (RUN/d).mkdir(parents=True,exist_ok=True)
def sha_bytes(b):return hashlib.sha256(b).hexdigest()
def snapshot(src,dest):
 if dest.exists():assert src.read_bytes()==dest.read_bytes(),f'Existing source snapshot differs: {dest}'
 else:shutil.copy2(src,dest)
 return {'source_file':str(src.relative_to(ROOT)),'snapshot_file':str(dest.relative_to(RUN)),'sha256':sha_bytes(dest.read_bytes())}
NUM=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
RESULT=re.compile(r'^\s*REMARK\s+VINA\s+RESULT\s*:\s*('+NUM+r')\s+('+NUM+r')\s+('+NUM+r')\s*$')
MARKER=re.compile(r'^\s*REMARK\s+VINA\s+RESULT\b')
MODEL=re.compile(r'^\s*MODEL\s+(\d+)\s*$')
def parse_pdbqt(text,drug,receptor,seed,source_file):
 records=[];issues=[];rank=None;model_ids=[];seen_in_model=0;model_line=None
 def issue(code,line,text):issues.append({'drug':drug,'receptor_block':receptor,'seed':seed,'source_file':source_file,'line_number':line,'issue':code,'detail':text})
 for line_no,line in enumerate(text.splitlines(),1):
  if re.match(r'^\s*MODEL\b',line):
   if rank is not None and seen_in_model!=1:issue('model_result_count',model_line,f'MODEL {rank}: {seen_in_model} VINA result records')
   m=MODEL.fullmatch(line);rank=int(m[1]) if m else None;seen_in_model=0;model_line=line_no
   if rank is None or rank<1:issue('malformed_model_rank',line_no,line);rank=None
   else:model_ids.append(rank)
  elif re.match(r'^\s*ENDMDL\b',line):
   if rank is not None and seen_in_model!=1:issue('model_result_count',model_line,f'MODEL {rank}: {seen_in_model} VINA result records')
   rank=None;seen_in_model=0;model_line=None
  elif MARKER.match(line):
   seen_in_model+=1;m=RESULT.fullmatch(line);nums=[np.nan]*3;status='valid'
   if m:
    nums=list(map(float,m.groups()))
    if not all(math.isfinite(v) for v in nums) or nums[1]<0 or nums[2]<nums[1]:status='malformed_or_invalid_numeric_record'
   else:status='malformed_vina_result'
   if status!='valid':issue(status,line_no,line)
   if rank is None:issue('result_without_model_rank',line_no,line)
   records.append({'drug':drug,'receptor_block':receptor,'seed':seed,'pose_rank':rank,'vina_score_kcal_mol':nums[0],'rmsd_lb':nums[1],'rmsd_ub':nums[2],'source_file':source_file,'line_number':line_no,'record_index':len(records)+1,'record_status':status,'raw_record':line})
 if rank is not None and seen_in_model!=1:issue('model_result_count',model_line,f'MODEL {rank}: {seen_in_model} records at EOF')
 if not records:issue('no_vina_result_records',None,'No REMARK VINA RESULT found')
 if len(model_ids)!=len(set(model_ids)):issue('duplicate_model_rank',None,str(model_ids))
 if model_ids and model_ids!=list(range(1,len(model_ids)+1)):issue('noncontiguous_or_unsorted_model_ranks',None,str(model_ids))
 return records,issues,len(model_ids)

def main():
 provenance=[]
 for src,dst in [(SOURCE/'config/docking_config.json',RUN/'config/docking_config_source.json'),(SOURCE/'tables/docking_run_manifest.csv',RUN/'config/docking_run_manifest_source.csv'),(SOURCE/'tables/raw_archive_manifest.csv',RUN/'config/raw_archive_manifest_source.csv'),(ROOT/'ANALYSIS_CONTRACT.md',RUN/'config/ANALYSIS_CONTRACT_source.md'),(ORDER_SOURCE,RUN/'config/drug_display_order.csv')]:provenance.append(snapshot(src,dst))
 config=json.loads((RUN/'config/docking_config_source.json').read_text());seeds=config['seeds'];assert len(seeds)==5 and len(set(seeds))==5
 drugs=pd.read_csv(RUN/'config/drug_display_order.csv').drug.tolist();receptors=['alpha1','alpha2'];assert len(drugs)==10
 dm=pd.read_csv(RUN/'config/docking_run_manifest_source.csv').set_index(['drug_id','receptor','seed'])
 am=pd.read_csv(RUN/'config/raw_archive_manifest_source.csv').set_index(['archive','member'])
 runmeta=json.loads((SOURCE/'run_manifest.json').read_text());allrows=[];allissues=[];files=[];archives=[]
 for drug in drugs:
  archive=SOURCE/'raw/archives'/f'{drug}.tar.gz';relative=str(archive.relative_to(ROOT));archive_sha=sha_bytes(archive.read_bytes()) if archive.exists() else None
  expected_hash=runmeta['artifact_sha256'].get(f'raw/archives/{drug}.tar.gz')
  if not archive.exists():
   archives.append({'archive':relative,'status':'missing','sha256':None,'expected_sha256':expected_hash});continue
  archive_ok=archive_sha==expected_hash
  archives.append({'archive':relative,'status':'hash_match' if archive_ok else 'hash_mismatch','sha256':archive_sha,'expected_sha256':expected_hash})
  with tarfile.open(archive,'r:gz') as tf:
   for member in tf.getmembers():
    if not member.isfile() or not member.name.endswith('.pdbqt'):continue
    source_file=relative+'::'+member.name
    mm=re.fullmatch(r'([^/]+)/(alpha[12])/(\d+)/([^/]+\.pdbqt)',member.name)
    if not mm:
     allissues.append({'drug':drug,'receptor_block':None,'seed':None,'source_file':source_file,'line_number':None,'issue':'unmapped_pdbqt_member','detail':member.name});continue
    d,rec,seed,filename=mm.groups();seed=int(seed);payload=tf.extractfile(member).read();digest=sha_bytes(payload)
    dest=RUN/'data/pdbqt'/d/rec/str(seed)/filename;dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():assert dest.read_bytes()==payload
    else:dest.write_bytes(payload)
    rows,issues,nmodels=parse_pdbqt(payload.decode('utf-8'),d,rec,seed,source_file)
    key=(d,rec,seed);known=key in dm.index;runhash=dm.loc[key,'output_sha256'] if known else None
    amkey=(f'raw/archives/{drug}.tar.gz',member.name);memberhash=am.loc[amkey,'sha256'] if amkey in am.index else None
    for condition,label in [(not archive_ok,'archive_hash_mismatch'),(digest!=runhash,'docking_output_hash_mismatch'),(digest!=memberhash,'archive_member_hash_mismatch'),(d!=drug,'archive_drug_mismatch'),(not known,'unexpected_drug_receptor_seed')]:
     if condition:issues.append({'drug':d,'receptor_block':rec,'seed':seed,'source_file':source_file,'line_number':None,'issue':label,'detail':digest})
    expected_poses=int(dm.loc[key,'n_poses']) if known else None
    if expected_poses is not None and len(rows)!=expected_poses:issues.append({'drug':d,'receptor_block':rec,'seed':seed,'source_file':source_file,'line_number':None,'issue':'pose_count_mismatch','detail':f'parsed={len(rows)}, manifest={expected_poses}'})
    files.append({'drug':d,'receptor_block':rec,'seed':seed,'source_file':source_file,'snapshot_file':str(dest.relative_to(RUN)),'sha256':digest,'docking_manifest_sha256':runhash,'archive_member_sha256':memberhash,'n_models':nmodels,'n_result_records':len(rows),'n_valid_records':sum(r['record_status']=='valid' for r in rows),'expected_n_poses':expected_poses,'n_issues':len(issues),'file_status':'ok' if not issues else 'flagged'})
    allrows.extend(rows);allissues.extend(issues)
 raw=pd.DataFrame(allrows);filedf=pd.DataFrame(files)
 rankmap={d:i for i,d in enumerate(drugs)}
 raw['drug_order']=raw.drug.map(rankmap);raw=raw.sort_values(['drug_order','receptor_block','seed','source_file','record_index']).drop(columns='drug_order')
 raw['source_pdbqt']=raw.source_file;raw['source_run']=str(SOURCE.relative_to(ROOT))
 raw.to_csv(RUN/'data/vina_scores_all_poses.csv',index=False,na_rep='NA')
 seedrows=[];missing=[]
 for drug in drugs:
  for rec in receptors:
   for seed in seeds:
    fs=filedf[(filedf.drug==drug)&(filedf.receptor_block==rec)&(filedf.seed==seed)]
    rr=raw[(raw.drug==drug)&(raw.receptor_block==rec)&(raw.seed==seed)]
    row={'drug':drug,'receptor_block':rec,'seed':seed,'pose_rank':np.nan,'vina_score_kcal_mol':np.nan,'rmsd_lb':np.nan,'rmsd_ub':np.nan,'source_file':None,'rank1_score_kcal_mol':np.nan,'rank1_is_best':None,'n_best_ties':0,'n_poses_preserved':len(rr),'status':'missing_output'}
    if len(fs)==0:missing.append({'drug':drug,'receptor_block':rec,'seed':seed,'issue':'missing_output'})
    elif len(fs)>1:row['status']='ambiguous_multiple_outputs';allissues.append({'drug':drug,'receptor_block':rec,'seed':seed,'source_file':' | '.join(fs.source_file),'line_number':None,'issue':'ambiguous_multiple_outputs','detail':f'{len(fs)} files; no seed score selected'})
    elif fs.iloc[0].file_status!='ok':row['status']='flagged_file_not_summarized';row['source_file']=fs.iloc[0].source_file
    else:
     bestscore=rr.vina_score_kcal_mol.min();ties=rr[rr.vina_score_kcal_mol==bestscore].sort_values(['pose_rank','record_index']);best=ties.iloc[0];rank1=rr[rr.pose_rank==1]
     rank1score=rank1.iloc[0].vina_score_kcal_mol if len(rank1)==1 else np.nan
     row.update({c:best[c] for c in ['pose_rank','vina_score_kcal_mol','rmsd_lb','rmsd_ub','source_file']})
     row.update({'rank1_score_kcal_mol':rank1score,'rank1_is_best':bool(rank1score==bestscore),'n_best_ties':len(ties),'status':'ok' if rank1score==bestscore else 'best_score_differs_from_rank1'})
     if rank1score!=bestscore:allissues.append({'drug':drug,'receptor_block':rec,'seed':seed,'source_file':best.source_file,'line_number':int(best.line_number),'issue':'rank1_best_score_mismatch','detail':'Seed best uses minimum valid score, lowest pose rank breaks ties'})
    seedrows.append(row)
 sb=pd.DataFrame(seedrows);sb['best_pose_rank']=sb.pose_rank;sb['best_vina_score_kcal_mol']=sb.vina_score_kcal_mol;sb['source_pdbqt']=sb.source_file;sb.to_csv(RUN/'data/vina_scores_seed_best.csv',index=False,na_rep='NA')
 summary=[]
 for drug in drugs:
  for rec in receptors:
   subset=sb[(sb.drug==drug)&(sb.receptor_block==rec)].set_index('seed');v=subset.vina_score_kcal_mol.dropna();n=len(v)
   out={'drug':drug,'receptor_block':rec,'median_best_seed_score':v.median() if n else np.nan,'mean_best_seed_score':v.mean() if n else np.nan,'min_best_seed_score':v.min() if n else np.nan,'max_best_seed_score':v.max() if n else np.nan,'IQR_best_seed_score':v.quantile(.75,interpolation='linear')-v.quantile(.25,interpolation='linear') if n else np.nan,'n_seeds_available':n,'n_seeds_expected':len(seeds),'sd_best_seed_score':v.std(ddof=1) if n>=2 else np.nan,'vina_favorability':-v.median() if n else np.nan}
   for seed in seeds:out[f'best_score_seed_{seed}']=subset.loc[seed,'vina_score_kcal_mol']
   summary.append(out)
 summary=pd.DataFrame(summary);summary.to_csv(RUN/'data/vina_scores_drug_receptor_summary.csv',index=False,na_rep='NA')
 filedf.to_csv(RUN/'data/pdbqt_file_manifest.csv',index=False,na_rep='NA');pd.DataFrame(archives).to_csv(RUN/'data/archive_manifest.csv',index=False,na_rep='NA')
 pd.DataFrame(allissues,columns=['drug','receptor_block','seed','source_file','line_number','issue','detail']).to_csv(RUN/'data/parsing_issues.csv',index=False,na_rep='NA')
 pd.DataFrame(missing,columns=['drug','receptor_block','seed','issue']).to_csv(RUN/'data/missing_seeds.csv',index=False)
 (RUN/'config/input_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')
 (RUN/'config/recovery_rules.json').write_text(json.dumps({'expected_seeds':seeds,'best_score':'minimum Vina docking score within seed; lowest pose_rank resolves ties; rank1 independently checked','IQR':'75th percentile minus 25th percentile; linear interpolation (NumPy/pandas type 7)','malformed':'all recognizable result lines retained with raw text/status; malformed records or integrity problems flag entire file, excluded from seed summary','missing':'expected missing seed retained as NA; never zero; summaries use available valid seeds','receptors':'separate alpha1 and alpha2; no combined score','clinical_data':'only archived drug-order CSV read; no phenotype values or optimization'},indent=2)+'\n')
 (RUN/'VINA_LAYER_FROZEN.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':{f:sha_bytes((RUN/'data'/f).read_bytes()) for f in ['vina_scores_all_poses.csv','vina_scores_seed_best.csv','vina_scores_drug_receptor_summary.csv']},'primary_score':'median of seed-level minima, by receptor independently','favorability':'negative of median raw score','clinical_values_used':False},indent=2)+'\n')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
 vals=summary.pivot(index='drug',columns='receptor_block',values='median_best_seed_score').reindex(index=drugs,columns=receptors);counts=summary.pivot(index='drug',columns='receptor_block',values='n_seeds_available').reindex_like(vals)
 fig,ax=plt.subplots(figsize=(11,8));cm=plt.get_cmap('viridis_r').copy();cm.set_bad('#d3d3d3');ax.set_facecolor('#d3d3d3');im=ax.pcolormesh(np.ma.masked_invalid(vals.values),cmap=cm,edgecolors='white',linewidth=.7)
 ax.set_ylim(10,0);ax.set_xticks([.5,1.5],receptors);ax.set_yticks(np.arange(10)+.5,drugs)
 for i in range(10):
  for j in range(2):
   val=vals.iloc[i,j];label=f'{val:.3f}\nn={int(counts.iloc[i,j])} seeds' if pd.notna(val) else f'NA\nn={int(counts.iloc[i,j])} seeds';color='white' if pd.notna(val) and im.norm(val)>.55 else '#111111';ax.text(j+.5,i+.5,label,ha='center',va='center',fontsize=11,color=color)
 fig.colorbar(im,ax=ax,label='Median best-seed Vina docking score (kcal/mol)');ax.set_title('Archived Vina docking scores across α1 and α2 receptor blocks\nMedian of independent seed-level best scores; no redocking',fontsize=14,pad=16)
 fig.text(.03,.02,'Computational docking scores, not experimentally measured binding affinity.\nFive search seeds per condition; not biological replicates.',fontsize=10);fig.tight_layout(rect=[0,.065,1,1])
 for ext in ['png','svg']:fig.savefig(RUN/'figures'/f'vina_score_alpha1_alpha2_heatmap.{ext}',dpi=200)
 plt.close(fig)
 fig,axs=plt.subplots(2,1,figsize=(13,10),sharex=True)
 for ax,rec in zip(axs,receptors):
  for j,seed in enumerate(seeds):
   v=sb[(sb.receptor_block==rec)&(sb.seed==seed)].set_index('drug').reindex(drugs).vina_score_kcal_mol
   ax.scatter(np.arange(10)+(j-2)*.055,v,s=38,label=str(seed),alpha=.8,zorder=3)
  for i,drug in enumerate(drugs):
   med=summary[(summary.drug==drug)&(summary.receptor_block==rec)].median_best_seed_score.iloc[0]
   if pd.notna(med):ax.hlines(med,i-.18,i+.18,color='#252d34',lw=1.8,zorder=4)
  ax.set_title('α1 receptor block' if rec=='alpha1' else 'α2 receptor block',fontsize=14,loc='left')
  ax.set_ylabel('Vina docking score (kcal/mol)');ax.grid(axis='y',alpha=.2)
 axs[1].set_xticks(range(10),drugs,rotation=35,ha='right')
 fig.suptitle('Archived Vina docking-score seed robustness',fontsize=18,y=.985)
 fig.legend(*axs[0].get_legend_handles_labels(),title='Independent computational search seed; black bar = median',loc='upper center',bbox_to_anchor=(.5,.95),ncol=5,fontsize=10)
 fig.text(.06,.016,'All five seed-best scores shown per drug/receptor. Computational search repeats, not biological replicates.',fontsize=11)
 fig.tight_layout(rect=[0,.045,1,.86])
 for ext in ['png','svg']:fig.savefig(RUN/'figures'/f'vina_score_seed_robustness.{ext}',dpi=200)
 plt.close(fig)
 qc={'docking_rerun':False,'PLIF_recomputed':False,'clinical_optimization':False,'drug_order':drugs,'receptor_blocks':receptors,'expected_seeds':seeds,'expected_outputs':len(drugs)*len(receptors)*len(seeds),'parsed_pdbqt_outputs':len(files),'preserved_result_records':len(raw),'valid_result_records':int((raw.record_status=='valid').sum()),'seed_rows':len(sb),'usable_seed_scores':int(sb.vina_score_kcal_mol.notna().sum()),'drug_receptor_summary_rows':len(summary),'missing_seeds':len(missing),'parsing_or_provenance_issues':len(allissues),'rank1_best_agreement':int(sb.rank1_is_best.fillna(False).astype(bool).sum()),'all_archives_hash_match':all(a['status']=='hash_match' for a in archives),'all_output_files_hash_match':all(f['sha256']==f['docking_manifest_sha256']==f['archive_member_sha256'] for f in files)}
 (RUN/'QC_REPORT.json').write_text(json.dumps(qc,indent=2)+'\n');print(json.dumps(qc,indent=2))


def analyze_layer():
 # Independent score layer is saved/frozen before any clinical fingerprint is opened.
 frozen=json.loads((RUN/'VINA_LAYER_FROZEN.json').read_text())
 for name,digest in frozen['files'].items():assert sha_bytes((RUN/'data'/name).read_bytes())==digest
 bridge=ROOT/'runs/20260920_clinical_bridge_hypothesis_poc'
 lineage=[]
 for src,dest in [(bridge/'raw/lineage/structural_fingerprint.csv',RUN/'data/frozen_plif_fingerprint.csv'),(bridge/'tables/clinical_fingerprint_10drug.csv',RUN/'data/frozen_clinical_fingerprint.csv'),(bridge/'STRUCTURAL_RESULTS_FROZEN.json',RUN/'config/STRUCTURAL_RESULTS_FROZEN_source.json'),(bridge/'CLINICAL_RESULTS_FROZEN.json',RUN/'config/CLINICAL_RESULTS_FROZEN_source.json'),(ROOT/'runs/20260920_archived_vina_score_recovery/scripts/recover_vina_scores.py',RUN/'config/original_recovery_script.py')]:lineage.append(snapshot(src,dest))
 sf=json.loads((RUN/'config/STRUCTURAL_RESULTS_FROZEN_source.json').read_text());cf=json.loads((RUN/'config/CLINICAL_RESULTS_FROZEN_source.json').read_text())
 assert sha_bytes((RUN/'data/frozen_plif_fingerprint.csv').read_bytes())==sf['sha256']
 assert sha_bytes((RUN/'data/frozen_clinical_fingerprint.csv').read_bytes())==cf['clinical_sha256']
 (RUN/'config/frozen_input_manifest.json').write_text(json.dumps(lineage,indent=2)+'\n')
 config={'clinical_terms':['Ataxia','Balance disorder','Coordination abnormal','Gait disturbance'],'primary_score':'median_best_seed_score','derived_score':'vina_favorability = -median_best_seed_score','missingness':'pairwise complete drugs for each term; no imputation','correlation':'Spearman across drugs; no p-values or hit calling; n<5 marked exploratory-insufficient; constant vectors return NA','PLIF_integration':'identity/source index only; no vector weighting, joint metric or feature selection','sd_definition':'sample SD of seed best scores, ddof=1','IQR_definition':'Q75-Q25, linear interpolation','interpretation':'Vina docking-score proxy, not experimental affinity; association/correspondence only'}
 (RUN/'config/layer_analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
 # No rule depends on observed clinical correlations.
 y=pd.read_csv(RUN/'data/frozen_clinical_fingerprint.csv',index_col=0);x=pd.read_csv(RUN/'data/frozen_plif_fingerprint.csv',index_col=0)
 drugs=pd.read_csv(RUN/'config/drug_display_order.csv').drug.tolist();assert list(y.index)==drugs and set(x.index)==set(drugs)
 summary=pd.read_csv(RUN/'data/vina_scores_drug_receptor_summary.csv')
 rawwide=summary.pivot(index='drug',columns='receptor_block',values='median_best_seed_score').reindex(index=drugs,columns=['alpha1','alpha2'])
 def corr(a,b):
  mask=a.notna()&b.notna();n=int(mask.sum());rho=np.nan;status='estimable'
  if n<5:status='exploratory-insufficient'
  elif a[mask].nunique()<2 or b[mask].nunique()<2:status='constant-vector'
  else:rho=float(spearmanr(a[mask],b[mask]).statistic)
  return rho,n,status,mask
 crossrho,crossn,crossstatus,_=corr(rawwide.alpha1,rawwide.alpha2)
 pd.DataFrame([{'x':'alpha1 median Vina score','y':'alpha2 median Vina score','n_effective':crossn,'spearman_rho':crossrho,'status':crossstatus}]).to_csv(RUN/'data/vina_alpha1_alpha2_correspondence.csv',index=False,na_rep='NA')
 output=[];joined=[]
 for rec in ['alpha1','alpha2']:
  for term in config['clinical_terms']:
   raw=rawwide[rec];clinical=y.loc[drugs,term];rr,n,status,mask=corr(raw,clinical);rf,nf,statusf,_=corr(-raw,clinical)
   assert n==nf and status==statusf
   if np.isfinite(rr):assert np.isclose(rr,-rf,atol=1e-14)
   output.append({'receptor_block':rec,'clinical_term':term,'n_effective':n,'rho_raw_vina_score':rr,'rho_vina_favorability':rf,'missing_drug_count':len(drugs)-n,'included_drugs':json.dumps(list(raw.index[mask])),'excluded_drugs':json.dumps(list(raw.index[~mask])),'status':status})
   for d in drugs:joined.append({'drug':d,'receptor_block':rec,'clinical_term':term,'median_best_seed_vina_score':raw.loc[d],'vina_favorability':-raw.loc[d],'clinical_logROR':clinical.loc[d],'included_in_correlation':bool(mask.loc[d])})
 out=pd.DataFrame(output);out.to_csv(RUN/'data/vina_clinical_c1_correspondence.csv',index=False,na_rep='NA')
 pd.DataFrame(joined).to_csv(RUN/'data/vina_clinical_c1_drug_values.csv',index=False,na_rep='NA')
 # PLIF references are exact identifiers, never new numerical summaries.
 descriptors={};index=[]
 plif_path='data/frozen_plif_fingerprint.csv';plifhash=sf['sha256']
 for d in drugs:
  row={'drug':d}
  for rec in ['alpha1','alpha2']:
   vectorid=f'sha256:{plifhash}#drug={d}&receptor_block={rec}'
   features=[f for f in x.columns if f.startswith(rec+'|')]
   descriptors[vectorid]={'source_file':plif_path,'source_sha256':plifhash,'row_key_column':'drug_id','row_key':d,'receptor_block':rec,'feature_names_in_frozen_order':features,'n_features':len(features)}
   row[f'{rec}_plif_vector_identifier']=vectorid
   row[f'{rec}_plif_source_file']=plif_path
   row[f'{rec}_plif_n_features']=len(features)
   row[f'{rec}_median_vina_score']=rawwide.loc[d,rec]
   row[f'{rec}_vina_favorability']=-rawwide.loc[d,rec]
  index.append(row)
 pd.DataFrame(index).to_csv(RUN/'data/structural_multiview_index.csv',index=False,na_rep='NA')
 (RUN/'data/plif_vector_reference_manifest.json').write_text(json.dumps(descriptors,indent=2)+'\n')
 plt.rcParams.update({'svg.fonttype':'none','font.family':'DejaVu Sans'})
 fig,ax=plt.subplots(figsize=(9,7))
 ax.scatter(rawwide.alpha1,rawwide.alpha2,color='#167d91',s=65,zorder=3)
 offsets={'diazepam':(-10,-18),'alprazolam':(8,-16),'triazolam':(8,12),'zolpidem':(8,10),'lorazepam':(-8,-16),'clonazepam':(8,8),'midazolam':(-8,14),'temazepam':(-8,12),'zopiclone':(8,-16),'zaleplon':(8,10)}
 for d in drugs:
  off=offsets[d];ax.annotate(d,(rawwide.loc[d,'alpha1'],rawwide.loc[d,'alpha2']),xytext=off,textcoords='offset points',ha='right' if off[0]<0 else 'left',fontsize=11)
 ax.set(xlabel='α1 median best-seed Vina docking score (kcal/mol)',ylabel='α2 median best-seed Vina docking score (kcal/mol)')
 ax.set_title(f'Archived Vina scores: α1 versus α2\nDescriptive Spearman rho = {crossrho:.3f}; n = {crossn} drugs',fontsize=15,pad=16)
 ax.margins(.16,.14);ax.grid(alpha=.15)
 fig.text(.04,.018,'More negative = more favorable docking score. Not receptor equivalence or experimental affinity.',fontsize=10)
 fig.tight_layout(rect=[0,.05,1,1])
 for ext in ['png','svg']:fig.savefig(RUN/'figures'/f'vina_score_alpha1_vs_alpha2.{ext}',dpi=200)
 plt.close(fig)
 cr=out.pivot(index='receptor_block',columns='clinical_term',values='rho_vina_favorability').reindex(index=['alpha1','alpha2'],columns=config['clinical_terms'])
 cn=out.pivot(index='receptor_block',columns='clinical_term',values='n_effective').reindex_like(cr)
 fig,ax=plt.subplots(figsize=(11,5.6));cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#d3d3d3')
 im=ax.pcolormesh(np.ma.masked_invalid(cr.values),cmap=cm,vmin=-1,vmax=1,edgecolors='white',linewidth=1)
 ax.set_ylim(2,0);ax.set_xticks(np.arange(4)+.5,config['clinical_terms'],rotation=15,ha='right');ax.set_yticks([.5,1.5],['α1','α2'])
 for i in range(2):
  for j in range(4):ax.text(j+.5,i+.5,f'n = {int(cn.iloc[i,j])}',ha='center',va='center',fontsize=14,color='white' if abs(cr.iloc[i,j])>.65 else '#222222')
 fig.colorbar(im,ax=ax,label='Spearman rho: Vina favorability vs clinical logROR',ticks=[-1,-.5,0,.5,1])
 ax.set_title('Vina docking-score correspondence with C1 Balance / Ataxia phenotype',fontsize=14,pad=34)
 ax.text(.5,1.035,'Favorability = −median Vina score; each cell compares drugs',ha='center',transform=ax.transAxes,fontsize=11)
 fig.text(.04,.02,'Observational unit = drug.\nVina score is a docking-score proxy, not experimental affinity.\nAssociation/correspondence only.',fontsize=11,linespacing=1.35)
 fig.tight_layout(rect=[0,.18,1,1])
 for ext in ['png','svg']:fig.savefig(RUN/'figures'/f'vina_clinical_c1_correspondence.{ext}',dpi=200)
 plt.close(fig)
 qc=json.loads((RUN/'QC_REPORT.json').read_text());qc.update({'frozen_X_sha256_verified':sf['sha256'],'frozen_Y_sha256_verified':cf['clinical_sha256'],'score_layer_frozen_before_clinical_correspondence':True,'sample_SD_ddof':1,'cross_receptor_descriptive_rho':crossrho,'C1_correspondence_rows':len(out),'raw_favorability_rho_sign_reversal_verified':True,'multiview_index_rows':len(index),'PLIF_vector_reference_count':len(descriptors),'PLIF_numeric_integration':False,'p_value_hit_calling':False})
 (RUN/'QC_REPORT.json').write_text(json.dumps(qc,indent=2)+'\n')
 print('Cross receptor rho:',crossrho);print(out[['receptor_block','clinical_term','n_effective','rho_raw_vina_score','rho_vina_favorability']].to_string(index=False))
if __name__=='__main__':
 main()
 analyze_layer()
