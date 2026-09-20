#!/usr/bin/env python3
"""Render archived X, Y, correlations and rotarod values without recomputation."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUN=Path(__file__).resolve().parents[1]
ROOT=RUN.parents[1]
SOURCE=ROOT/'runs/20260920_clinical_bridge_hypothesis_poc'
DATA=RUN/'data';FIG=RUN/'figures'
for p in [DATA,FIG]:p.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files={
 'ANALYSIS_CONTRACT_source.md':'config/ANALYSIS_CONTRACT_source.md',
 'structural_fingerprint.csv':'raw/lineage/structural_fingerprint.csv',
 'clinical_fingerprint_10drug.csv':'tables/clinical_fingerprint_10drug.csv',
 'structural_feature_clinical_correlations.csv':'tables/structural_feature_clinical_correlations.csv',
 'invivo_anchor_contexts.csv':'tables/invivo_anchor_contexts.csv',
 'clinical_term_mapping.csv':'tables/clinical_term_mapping.csv',
 'STRUCTURAL_RESULTS_FROZEN.json':'STRUCTURAL_RESULTS_FROZEN.json',
 'CLINICAL_RESULTS_FROZEN.json':'CLINICAL_RESULTS_FROZEN.json',
 'structural_fingerprint_frozen_manifest.csv':'tables/structural_fingerprint_frozen_manifest.csv',
 'clinical_phenotype_config.json':'config/clinical_phenotype_config.json',
}
manifest=[]
for name,rel in files.items():
 src=SOURCE/rel;dest=DATA/name
 if src.exists():
  if dest.exists():assert sha(dest)==sha(src),f'Immutable input mismatch: {name}'
  else:shutil.copy2(src,dest)
  manifest.append({'file':name,'source':str(src.relative_to(ROOT)),'sha256':sha(dest)})
 elif not dest.exists():raise FileNotFoundError(src)
# Permit standalone regeneration from the bundled, unmodified data snapshots.
if manifest:(RUN/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
sf=json.loads((DATA/'STRUCTURAL_RESULTS_FROZEN.json').read_text());yf=json.loads((DATA/'CLINICAL_RESULTS_FROZEN.json').read_text())
assert sha(DATA/'structural_fingerprint.csv')==sf['sha256']
assert sha(DATA/'clinical_fingerprint_10drug.csv')==yf['clinical_sha256']
x=pd.read_csv(DATA/'structural_fingerprint.csv',index_col=0)
y=pd.read_csv(DATA/'clinical_fingerprint_10drug.csv',index_col=0)
assert x.index.is_unique and y.index.is_unique and set(x.index)==set(y.index)
# Display order is exactly the existing frozen Y order, never a clustering result.
order=list(y.index);x=x.loc[order];assert x.shape==(10,48) and y.shape==(10,14)
assert x.columns.str.startswith('alpha1|').sum()==21
corr=pd.read_csv(DATA/'structural_feature_clinical_correlations.csv')
assert not corr.duplicated(['structural_feature','phenotype_term']).any()
assert len(corr)==48*14
rho=corr.pivot(index='structural_feature',columns='phenotype_term',values='rho').reindex(index=x.columns,columns=y.columns)
ne=corr.pivot(index='structural_feature',columns='phenotype_term',values='n_effective').reindex_like(rho)
assert ne.notna().all().all()
anchor=pd.read_csv(DATA/'invivo_anchor_contexts.csv')
assert set(anchor.drug)=={'diazepam','triazolam'} and set(anchor.study_id)=={'Nishino2008'}
assert len(anchor)==16 and (anchor.tier==1).all()
for c in ['dose','observation_time','value']:anchor[c]=pd.to_numeric(anchor[c],errors='raise')
assert set(anchor.dose)=={2,5} and set(anchor.observation_time)=={15,30,60,90}
axes=json.loads((DATA/'clinical_phenotype_config.json').read_text())['axes']
term_labels=[next(a for a,tt in axes.items() if t in tt)+' '+t for t in y.columns]
feature_ids=[f'F{i+1:02}' for i in range(48)]
key=pd.DataFrame({'display_id':feature_ids,'frozen_feature':x.columns,'source_column_order':range(1,49)})
key.to_csv(DATA/'feature_display_key.csv',index=False)
pd.DataFrame({'display_row':range(1,11),'drug':order}).to_csv(DATA/'drug_display_order.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.facecolor':'white'})
fig=plt.figure(figsize=(24,26),facecolor='white')
gs=fig.add_gridspec(3,1,height_ratios=[5.5,13.2,3.5],left=.115,right=.95,top=.855,bottom=.055,hspace=.45)
fig.suptitle('Frozen structural and clinical fingerprints aligned by drug identity',x=.5,y=.983,fontsize=25,fontweight='bold')
fig.text(.5,.959,'association/correspondence only; not causation or mediation',ha='center',fontsize=16,color='#8c3e25')
# Panel A: equal row heights and independent color scales; no cross-feature arrows.
ga=gs[0].subgridspec(1,2,width_ratios=[1.55,1],wspace=.24)
ax=fig.add_subplot(ga[0]);ay=fig.add_subplot(ga[1])
xc=plt.get_cmap('viridis').copy();xc.set_bad('#bdbdbd')
yc=plt.get_cmap('RdBu_r').copy();yc.set_bad('#bdbdbd')
imx=ax.imshow(x.values,aspect='auto',cmap=xc,vmin=0,vmax=1,interpolation='nearest')
# Untransformed frozen values; symmetric color range encompasses all finite Y values.
ylim=float(np.ceil(np.nanmax(np.abs(y.values))))
imy=ay.imshow(y.values,aspect='auto',cmap=yc,vmin=-ylim,vmax=ylim,interpolation='nearest')
for a in [ax,ay]:
 a.set_yticks(range(10),order,fontsize=12);a.set_ylim(9.5,-.5)
 a.set_yticks(np.arange(-.5,10,1),minor=True);a.grid(which='minor',axis='y',color='white',lw=.45);a.tick_params(which='minor',length=0)
ax.set_xticks(range(48),feature_ids,rotation=90,fontsize=9)
ay.set_xticks(range(14),term_labels,rotation=55,ha='right',fontsize=10)
ax.axvline(20.5,color='white',lw=2)
fig.text(ax.get_position().x0,.886,'A  |  X: 48-feature structural PLIF',fontsize=17,fontweight='bold')
fig.text(ay.get_position().x0,.886,'Y: 14-term clinical fingerprint',fontsize=17,fontweight='bold')
ax.text(10/47,1.015,'alpha1 | F01-F21',ha='center',transform=ax.transAxes,fontsize=11)
ax.text(34/47,1.015,'alpha2 | F22-F48',ha='center',transform=ax.transAxes,fontsize=11)
for a,im,label in [(ax,imx,'Frozen contact frequency'),(ay,imy,'Frozen log(ROR) | all-role FAERS')]:
 box=a.get_position();cax=fig.add_axes([box.x0,.919,box.width,.006]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.ax.xaxis.set_ticks_position('top');cb.ax.xaxis.set_label_position('top');cb.set_label(label,fontsize=10,labelpad=5);cb.ax.tick_params(labelsize=9)
# Panel B reuses saved Spearman statistics; does not calculate correlations.
gb=gs[1].subgridspec(1,2,width_ratios=[1.9,.65],wspace=.13)
b=fig.add_subplot(gb[0]);notes=fig.add_subplot(gb[1]);notes.axis('off')
cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#bdbdbd')
im=b.imshow(rho.values,aspect='auto',vmin=-1,vmax=1,cmap=cm,interpolation='nearest')
flabels=[fid+'  '+f.replace('|any_contact_frequency','').replace('BZD_GAMMA2_','gamma2:').replace('BZD_SITE_','site:').replace('|', ' / ') for fid,f in zip(feature_ids,x.columns)]
b.set_yticks(range(48),flabels,fontsize=9)
b.set_xticks(range(14),term_labels,rotation=50,ha='right',fontsize=10)
b.axhline(20.5,color='#343434',lw=1.5)
b.set_title('B  |  Archived drug-level Spearman correlation matrix',loc='left',fontsize=17,pad=18,fontweight='bold')
for i in range(48):
 for j in range(14):
  value=rho.iloc[i,j];color='white' if np.isfinite(value) and abs(value)>.68 else '#202020'
  b.text(j,i,str(int(ne.iloc[i,j])),ha='center',va='center',fontsize=8,color=color)
cax=notes.inset_axes([.10,.85,.8,.025]);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,-.5,0,.5,1]);cb.set_label('Spearman rho',fontsize=12)
notes.text(.08,.98,'HOW TO READ',fontsize=14,fontweight='bold',va='top')
notes.text(.08,.78,'Cell color: archived rho\nCell number: n_effective\nUnit of analysis: drug (n <= 10)\n\nAll 48 X features and 14 Y terms\nare retained in frozen order.\nNo correlation recalculation.\n\nGray = NA, never zero.\nY: event count <5 is masked\nin the existing frozen input.\nRho: undefined values remain NA.\n\nTop panels share drug identity\nand exactly the same row order.\nThey use separate color scales.\nNo joint clustering or rescaling.\n\nFAERS log(ROR) is a reporting\nsignal, not incidence or risk ratio.\nNo causal feature-term arrows.\n\nStructural caveat:\nalpha2 is a provisional construct;\nreceptor blocks are not replicates.',va='top',fontsize=12,linespacing=1.65)
# Panel C: raw archived fractions, separate dose panels; no two-drug correlation.
gc=gs[2].subgridspec(1,3,width_ratios=[1,1,.95],wspace=.3)
for k,dose in enumerate([2,5]):
 a=fig.add_subplot(gc[k])
 for drug,color in [('diazepam','#15798a'),('triazolam','#cc7935')]:
  sub=anchor[(anchor.drug==drug)&(anchor.dose==dose)].sort_values('observation_time')
  a.plot(sub.observation_time,sub.value,marker='o',lw=2,ms=7,color=color,label=drug)
 a.set(title=f'{dose} mg/kg, oral',xlabel='Minutes after administration',ylabel='Rotarod failure fraction',ylim=(-.04,1.04),xticks=[15,30,60,90]);a.grid(axis='y',alpha=.18);a.legend(frameon=False,loc='upper right',fontsize=10)
 if k==0:a.text(0,1.23,'C  |  External biological anchor: Nishino 2008',transform=a.transAxes,fontsize=17,fontweight='bold')
a=fig.add_subplot(gc[2]);a.axis('off');a.text(0,.95,'Diazepam / triazolam only\nMouse, ICR, male; n=10 mice per\nreported dose/time condition.\n\nArchived observed fractions;\nno dose pooling or normalization.\nRepeated dose/time conditions\nare not independent drug pairs.\n\nExternal anchor, not a mediator.\nNo n=2 correlation or causal path.\nDOI: 10.1254/jphs.08107FP',va='top',fontsize=12,linespacing=1.5)
fig.text(.115,.018,'Frozen X/Y are unchanged. Linkage key: generic drug identity. F01-F48 map to exact source features in data/feature_display_key.csv.',fontsize=12)
for ext in ['png','svg']:fig.savefig(FIG/f'frozen_xy_paired_heatmap_with_external_anchor.{ext}',dpi=180)
plt.close(fig)
# Report only representation validation, never modify numerical inputs.
qc={'X_frozen_sha256':sha(DATA/'structural_fingerprint.csv'),'Y_frozen_sha256':sha(DATA/'clinical_fingerprint_10drug.csv'),'X_shape':list(x.shape),'Y_shape':list(y.shape),'drug_order':order,'Y_NA_cells':int(y.isna().sum().sum()),'correlation_cells':int(rho.size),'n_effective_annotations':int(ne.size),'anchor_rows':len(anchor),'feature_selection':False,'input_recalculation':False,'correlation_recalculation':False,'same_drug_order':True,'causal_arrows':False,'source_run':str(SOURCE.relative_to(ROOT))}
(RUN/'QC_REPORT.json').write_text(json.dumps(qc,indent=2)+'\n')
print(json.dumps(qc,indent=2))
