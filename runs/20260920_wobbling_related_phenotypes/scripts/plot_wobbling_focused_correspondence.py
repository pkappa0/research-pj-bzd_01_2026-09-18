#!/usr/bin/env python3
"""Render frozen inputs and archived correlations; summarize C1 directions only."""
from pathlib import Path
import hashlib, json, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

RUN=Path(__file__).resolve().parents[1]
ROOT=RUN.parents[1]
SRC=ROOT/'runs/20260920_clinical_bridge_hypothesis_poc'
DATA=RUN/'data'; TABLES=RUN/'tables'; FIGURES=RUN/'figures'
for p in [DATA,TABLES,FIGURES]:p.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inputs={
 'structural_fingerprint.csv':'raw/lineage/structural_fingerprint.csv',
 'clinical_fingerprint_10drug.csv':'tables/clinical_fingerprint_10drug.csv',
 'structural_feature_clinical_correlations.csv':'tables/structural_feature_clinical_correlations.csv',
 'invivo_anchor_contexts.csv':'tables/invivo_anchor_contexts.csv',
 'STRUCTURAL_RESULTS_FROZEN.json':'STRUCTURAL_RESULTS_FROZEN.json',
 'CLINICAL_RESULTS_FROZEN.json':'CLINICAL_RESULTS_FROZEN.json',
 'clinical_phenotype_config.json':'config/clinical_phenotype_config.json',
 'ANALYSIS_CONTRACT_source.md':'config/ANALYSIS_CONTRACT_source.md',
}
old_manifest=RUN/'input_manifest.json'
if old_manifest.exists():
 for item in json.loads(old_manifest.read_text()):
  assert sha(DATA/item['file'])==item['sha256'],f"Input altered: {item['file']}"
manifest=[]
for name,rel in inputs.items():
 src=SRC/rel;dest=DATA/name
 if src.exists():
  if dest.exists():assert sha(src)==sha(dest),f'Source mismatch: {name}'
  else:shutil.copy2(src,dest)
 elif not dest.exists():raise FileNotFoundError(src)
 manifest.append({'file':name,'source':str((SRC/rel).relative_to(ROOT)),'sha256':sha(dest)})
old_manifest.write_text(json.dumps(manifest,indent=2)+'\n')
sf=json.loads((DATA/'STRUCTURAL_RESULTS_FROZEN.json').read_text());yf=json.loads((DATA/'CLINICAL_RESULTS_FROZEN.json').read_text())
assert sha(DATA/'structural_fingerprint.csv')==sf['sha256']
assert sha(DATA/'clinical_fingerprint_10drug.csv')==yf['clinical_sha256']
x=pd.read_csv(DATA/'structural_fingerprint.csv',index_col=0);y=pd.read_csv(DATA/'clinical_fingerprint_10drug.csv',index_col=0)
assert x.shape==(10,48) and y.shape==(10,14)
assert x.index.is_unique and y.index.is_unique and set(x.index)==set(y.index)
DRUGS=list(y.index); FEATURES=list(x.columns); IDS=[f'F{i+1:02}' for i in range(48)]; x=x.loc[DRUGS]
C1=['Ataxia','Balance disorder','Coordination abnormal','Gait disturbance']
C5=['Muscular weakness','Hypotonia']
CONTEXT=['Somnolence','Sedation','Dizziness','Vertigo','Fall']
AXIS={**dict.fromkeys(C1,'C1'),**dict.fromkeys(C5,'C5'),**dict.fromkeys(CONTEXT[:2],'C3'),**dict.fromkeys(CONTEXT[2:4],'C4'),'Fall':'C2'}
EXAMPLES=['F02','F04','F16','F28']
arch=pd.read_csv(DATA/'structural_feature_clinical_correlations.csv')
assert len(arch)==672 and not arch.duplicated(['structural_feature','phenotype_term']).any()
lookup=arch.set_index(['structural_feature','phenotype_term'])
def extract(terms):
 rows=[]
 for fid,f in zip(IDS,FEATURES):
  for t in terms:
   r=lookup.loc[(f,t)];assert r.phenotype_axis==AXIS[t]
   rows.append({'feature_id':fid,'feature_label':f,'clinical_axis':AXIS[t],'clinical_term':t,'spearman_rho':r.rho,'n_effective':int(r.n_effective)})
 return pd.DataFrame(rows)
core=extract(C1+C5);context=extract(CONTEXT)
core.to_csv(TABLES/'wobbling_plif_clinical_correspondence.csv',index=False,na_rep='NA')
context.to_csv(TABLES/'context_plif_clinical_correspondence.csv',index=False,na_rep='NA')
# Only new arithmetic: explicitly requested descriptive summaries of archived C1 rho.
summary=[]
for fid,f in zip(IDS,FEATURES):
 r=core[(core.feature_id==fid)&(core.clinical_axis=='C1')].spearman_rho.dropna()
 summary.append({'feature_id':fid,'feature_label':f,'number_estimable_c1_terms':len(r),'number_positive_rho':int((r>0).sum()),'number_negative_rho':int((r<0).sum()),'number_zero_rho':int((r==0).sum()),'median_rho_across_c1_terms':r.median() if len(r) else np.nan,'min_rho':r.min() if len(r) else np.nan,'max_rho':r.max() if len(r) else np.nan})
pd.DataFrame(summary).to_csv(TABLES/'c1_direction_consistency_summary.csv',index=False,na_rep='NA')
pd.DataFrame({'feature_id':IDS,'feature_label':FEATURES,'source_column_order':range(1,49)}).to_csv(TABLES/'feature_display_key.csv',index=False)
pd.DataFrame({'display_row':range(1,11),'drug':DRUGS}).to_csv(TABLES/'drug_display_order.csv',index=False)
example=x.loc[:,[FEATURES[IDS.index(fid)] for fid in EXAMPLES]].copy();example.columns=EXAMPLES
example.to_csv(TABLES/'representative_frozen_plif_values.csv',na_rep='NA')
anchor=pd.read_csv(DATA/'invivo_anchor_contexts.csv')
assert len(anchor)==16 and set(anchor.drug)=={'diazepam','triazolam'}
assert set(anchor.study_id)=={'Nishino2008'} and (anchor.tier==1).all()
assert set(anchor.dose)=={2,5} and set(anchor.observation_time)=={15,30,60,90}
assert not anchor.duplicated(['drug','dose','observation_time']).any()

def matrices(df,terms):
 return tuple(df.pivot(index='feature_id',columns='clinical_term',values=v).reindex(index=IDS,columns=terms) for v in ['spearman_rho','n_effective'])
rho,ncore=matrices(core,C1+C5);ctx,nctx=matrices(context,CONTEXT)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(24,28),facecolor='white')
gs=fig.add_gridspec(3,1,height_ratios=[5.3,13.6,4.4],left=.16,right=.94,top=.925,bottom=.115,hspace=.44)
fig.suptitle('Wobbling-focused PLIF-clinical phenotype correspondence',fontsize=26,fontweight='bold',y=.981)
fig.text(.5,.959,'association/correspondence only; not causation or mediation',ha='center',fontsize=16,color='#8b432c')
# A: categorical organization only. No arrows, path models or implied estimates.
a=fig.add_subplot(gs[0]);a.axis('off')
a.text(0,1.035,'A  |  Clinical phenotype definition',fontsize=18,fontweight='bold',transform=a.transAxes)
a.text(.5,.955,'WOBBLING / MOTOR INSTABILITY',ha='center',fontsize=19,fontweight='bold')
a.text(.5,.865,'Phenotype organization for exploratory correspondence analysis',ha='center',fontsize=12,color='#4f5963')
def box(ax,xx,yy,ww,hh,color):
 ax.add_patch(FancyBboxPatch((xx,yy),ww,hh,boxstyle='round,pad=0.01,rounding_size=0.015',facecolor=color,edgecolor='#bcc7ce',lw=1,transform=ax.transAxes,clip_on=False))
box(a,.01,.40,.48,.37,'#e7f2f5');box(a,.52,.40,.47,.37,'#f7eddf')
a.text(.035,.715,'PRIMARY PHENOTYPE  |  C1 Balance / Ataxia',fontsize=17,fontweight='bold')
a.text(.035,.635,'Ataxia\nBalance disorder\nCoordination abnormal\nGait disturbance',va='top',fontsize=13,linespacing=1.3)
a.text(.545,.715,'RELATED MOTOR PHENOTYPE  |  C5 Motor weakness',fontsize=15,fontweight='bold')
a.text(.545,.615,'Muscular weakness\nHypotonia',va='top',fontsize=14,linespacing=1.7)
a.text(.5,.325,'RELATED CLINICAL PHENOTYPES',ha='center',fontsize=14,fontweight='bold',color='#56616b')
for xx,title in [(.01,'C3 Sedation'),(.35,'C4 Dizziness'),(.69,'C2 Fall')]:
 box(a,xx,.025,.30,.22,'#f3f4f6');a.text(xx+.15,.125,title,ha='center',fontsize=16,fontweight='bold')
# B and E share feature order and dimensions, but remain separate clinical groups.
gmid=gs[1].subgridspec(1,2,width_ratios=[1.25,1],wspace=.23)
b=fig.add_subplot(gmid[0]);e=fig.add_subplot(gmid[1])
def matrix_plot(ax,values,nvalues,labels,rows):
 ax.set_facecolor('#bdbdbd');cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#bdbdbd')
 im=ax.pcolormesh(np.ma.masked_invalid(values.to_numpy(float)),cmap=cm,vmin=-1,vmax=1,edgecolors='none',rasterized=False)
 nr,nc=values.shape;ax.set_xlim(0,nc);ax.set_ylim(nr,0)
 ax.set_xticks(np.arange(nc)+.5,labels,rotation=38,ha='right',fontsize=11)
 ax.set_yticks(np.arange(nr)+.5,rows,fontsize=9.5);ax.tick_params(length=0)
 ax.axhline(21,color='#3b4046',lw=1.7)
 for i in range(nr):
  for j in range(nc):
   v=values.iloc[i,j];color='white' if pd.notna(v) and abs(v)>.68 else '#222222'
   ax.text(j+.5,i+.5,str(int(nvalues.iloc[i,j])),ha='center',va='center',fontsize=9,color=color)
 return im
labels=[fid+'  '+f.replace('|any_contact_frequency','').replace('BZD_GAMMA2_','gamma2:').replace('BZD_SITE_','site:').replace('|',' / ') for fid,f in zip(IDS,FEATURES)]
im=matrix_plot(b,rho,ncore,[AXIS[t]+' '+t for t in C1+C5],labels)
b.axvline(4,color='white',lw=7);b.axvline(4,color='#6a727a',lw=1.5)
b.set_title('B  |  Drug-level PLIF correspondence\nwith balance / motor phenotypes',loc='left',fontsize=18,fontweight='bold',pad=40)
b.text(2/6,1.013,'PRIMARY: C1',transform=b.transAxes,ha='center',fontsize=13,fontweight='bold',color='#206576')
b.text(5/6,1.013,'RELATED MOTOR: C5',transform=b.transAxes,ha='center',fontsize=13,fontweight='bold',color='#845a21')
matrix_plot(e,ctx,nctx,[AXIS[t]+' '+t for t in CONTEXT],IDS)
e.set_title('E  |  Related clinical phenotypes\nnot primary wobbling endpoint',loc='left',fontsize=17,fontweight='bold',pad=40)
e.text(.5,1.013,'RELATED CLINICAL: C3 / C4 / C2',transform=e.transAxes,ha='center',fontsize=13,color='#56616b')
# Shared legend explicitly describes both archived correlation matrices.
boxpos=e.get_position();cax=fig.add_axes([.955,boxpos.y0+.12, .012,boxpos.height*.55]);cb=fig.colorbar(im,cax=cax,ticks=[-1,-.5,0,.5,1]);cb.set_label('Archived Spearman rho | cell number = n_effective',fontsize=12)
# C / D lower panels preserve original observations. C examples are user-specified.
gbottom=gs[2].subgridspec(1,3,width_ratios=[.85,1,1],wspace=.40)
c=fig.add_subplot(gbottom[0]);cm=plt.get_cmap('viridis').copy();cm.set_bad('#bdbdbd');c.set_facecolor('#bdbdbd')
ic=c.pcolormesh(np.ma.masked_invalid(example.values),cmap=cm,vmin=0,vmax=1,rasterized=False)
c.set_ylim(10,0);c.set_xticks(np.arange(4)+.5,EXAMPLES,fontsize=12);c.set_yticks(np.arange(10)+.5,DRUGS,fontsize=11);c.tick_params(length=0)
c.set_title('C  |  Representative visually variable\nPLIF features',loc='left',fontsize=16,fontweight='bold',pad=30)
cbax=c.inset_axes([0,-.16,1,.035]);cb=fig.colorbar(ic,cax=cbax,orientation='horizontal',ticks=[0,.5,1]);cb.set_label('Frozen contact frequency',fontsize=11);cb.ax.tick_params(labelsize=10)
for k,dose in enumerate([2,5]):
 d=fig.add_subplot(gbottom[k+1])
 for drug,col in [('diazepam','#197b8b'),('triazolam','#cc7935')]:
  rows=anchor[(anchor.drug==drug)&(anchor.dose==dose)].sort_values('observation_time')
  d.plot(rows.observation_time,rows.value,marker='o',lw=2,ms=6,color=col,label=drug)
 d.set_title(f'{dose} mg/kg, oral',fontsize=14,pad=12)
 d.set(ylim=(-.04,1.04),xticks=[15,30,60,90],xlabel='Minutes after administration')
 if k==0:
  d.set_ylabel('Rotarod failure fraction');p=d.get_position();fig.text(p.x0,p.y1+.031,'D  |  External biological anchor: Nishino 2008',fontsize=17,fontweight='bold')
 d.grid(axis='y',alpha=.2);d.legend(loc='upper right',frameon=False,fontsize=11)
fig.text(.16,.043,'Panel C: pre-specified for visualization only; no feature selection.\nAll 48 features remain in B/E; examples are not biomarkers or hits.',fontsize=11,linespacing=1.6)
fig.text(.47,.043,'Context-matched diazepam/triazolam comparison only (n=2 drugs).\nExternal anchor, not a mediator. No correlation or causal inference.',fontsize=11,linespacing=1.6)
fig.text(.16,.015,'Each correlation is calculated across drugs. Rows are structural PLIF features, not causal residues. Gray = NA / undefined; not zero.',fontsize=12)
for ext in ['png','svg']:fig.savefig(FIGURES/f'wobbling_focused_plif_clinical_correspondence.{ext}',dpi=200)
plt.close(fig)
# Verify exported values against archived cells, without recomputing any correlation.
for frame in [core,context]:
 for r in frame.itertuples():
  original=lookup.loc[(r.feature_label,r.clinical_term)]
  assert (pd.isna(r.spearman_rho) and pd.isna(original.rho)) or r.spearman_rho==original.rho
  assert r.n_effective==original.n_effective
assert len(core)==288 and len(context)==240 and len(summary)==48
assert list(pd.DataFrame(summary).feature_id)==IDS
assert all(v['number_positive_rho']+v['number_negative_rho']+v['number_zero_rho']==v['number_estimable_c1_terms'] for v in summary)
qc={'X_sha256':sha(DATA/'structural_fingerprint.csv'),'Y_sha256':sha(DATA/'clinical_fingerprint_10drug.csv'),'source_hashes_verified':True,'frozen_X_shape':list(x.shape),'frozen_Y_shape':list(y.shape),'core_cells':len(core),'context_cells':len(context),'C1_summary_rows':len(summary),'n_effective_annotated_cells':len(core)+len(context),'feature_order':IDS,'drug_order':DRUGS,'core_term_order':C1+C5,'context_term_order':CONTEXT,'representative_features_user_specified':EXAMPLES,'Y_NA_cells_unchanged':int(y.isna().sum().sum()),'core_undefined_rho_cells':int(core.spearman_rho.isna().sum()),'anchor_rows':len(anchor),'all_exported_correlations_equal_archive':True,'new_calculations':'Only descriptive C1 sign counts, median, minimum and maximum of archived rho','PLIF_ROR_or_correlation_recalculation':False,'statistical_feature_selection_or_ranking':False,'causal_arrows':False}
(RUN/'QC_REPORT.json').write_text(json.dumps(qc,indent=2)+'\n')
print('Rendered figure; 288 core / 240 context cells; 48 unsorted C1 descriptive summaries. Source hashes unchanged.')
