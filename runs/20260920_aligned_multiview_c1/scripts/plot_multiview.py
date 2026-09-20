#!/usr/bin/env python3
"""Plot frozen views; no docking, PLIF, ROR, or correlation calculation."""
from pathlib import Path
import json, hashlib, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
R=Path(__file__).resolve().parents[1]; ROOT=R.parents[1]
D=R/'data'; F=R/'figures'
D.mkdir(exist_ok=True); F.mkdir(exist_ok=True)
SOURCES={
 'plif.csv':'runs/20260920_vina_score_layer/data/frozen_plif_fingerprint.csv',
 'clinical.csv':'runs/20260920_vina_score_layer/data/frozen_clinical_fingerprint.csv',
 'vina_summary.csv':'runs/20260920_vina_score_layer/data/vina_scores_drug_receptor_summary.csv',
 'plif_c1_correspondence.csv':'runs/20260920_wobbling_focused_correspondence/tables/wobbling_plif_clinical_correspondence.csv',
 'vina_c1_correspondence.csv':'runs/20260920_vina_score_layer/data/vina_clinical_c1_correspondence.csv',
 'feature_display_key.csv':'runs/20260920_wobbling_focused_correspondence/tables/feature_display_key.csv',
 'drug_display_order.csv':'runs/20260920_vina_score_layer/config/drug_display_order.csv'}
manifest=[]
for dest,src in SOURCES.items():
 b=(ROOT/src).read_bytes(); p=D/dest
 if p.exists(): assert p.read_bytes()==b
 else:p.write_bytes(b)
 manifest.append(dict(source=src,snapshot='data/'+dest,sha256=hashlib.sha256(b).hexdigest()))
b=(ROOT/'ANALYSIS_CONTRACT.md').read_bytes();(R/'ANALYSIS_CONTRACT_source.md').write_bytes(b)
manifest.append(dict(source='ANALYSIS_CONTRACT.md',snapshot='ANALYSIS_CONTRACT_source.md',sha256=hashlib.sha256(b).hexdigest()))
(R/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
drugs=['diazepam','alprazolam','triazolam','zolpidem','lorazepam','clonazepam','midazolam','temazepam','zopiclone','zaleplon']
terms=['Ataxia','Balance disorder','Coordination abnormal','Gait disturbance']; ids=[f'F{i:02}' for i in range(1,49)]
assert pd.read_csv(D/'drug_display_order.csv').drug.tolist()==drugs
key=pd.read_csv(D/'feature_display_key.csv');assert key.feature_id.tolist()==ids
x=pd.read_csv(D/'plif.csv').set_index('drug_id').loc[drugs]
assert x.columns.tolist()==key.feature_label.tolist()
assert all(c.startswith('alpha1|') for c in x.columns[:21]) and all(c.startswith('alpha2|') for c in x.columns[21:])
y=pd.read_csv(D/'clinical.csv').set_index('drug').loc[drugs,terms]
v=pd.read_csv(D/'vina_summary.csv').pivot(index='drug',columns='receptor_block',values='median_best_seed_score').loc[drugs,['alpha1','alpha2']]
# Only a display transform: rank original favorability within each receptor.
# Rank 1 = least favorable, rank 10 = most favorable; average ties.
ranks=(-v).rank(axis=0,method='average',ascending=True)
ranks.rename(columns=lambda c:c+'_within_receptor_favorability_rank').to_csv(D/'vina_display_ranks.csv')
pc=pd.read_csv(D/'plif_c1_correspondence.csv');vc=pd.read_csv(D/'vina_c1_correspondence.csv')
pr=pc.pivot(index='feature_id',columns='clinical_term',values='spearman_rho').loc[ids,terms]
pn=pc.pivot(index='feature_id',columns='clinical_term',values='n_effective').loc[ids,terms]
vr=vc.pivot(index='receptor_block',columns='clinical_term',values='rho_vina_favorability').loc[['alpha1','alpha2'],terms]
vn=vc.pivot(index='receptor_block',columns='clinical_term',values='n_effective').loc[['alpha1','alpha2'],terms]
assert pr.shape==(48,4) and vr.shape==(2,4)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(22,21),facecolor='white')
fig.suptitle('Multi-view structural correspondence with C1 Balance / Ataxia phenotype',fontsize=23,y=.985)
fig.text(.5,.955,'PLIF interaction pattern and Vina docking score are shown as\ncomplementary structural views aligned by drug identity',ha='center',va='top',fontsize=15)
fig.text(.07,.91,'A  Drug-level multi-view profile',fontsize=18,weight='bold')
def heat(ax,a,cmap,vmin,vmax):
 cm=plt.get_cmap(cmap).copy();cm.set_bad('#c5c5c5')
 im=ax.imshow(np.ma.masked_invalid(np.asarray(a,dtype=float)),aspect='auto',interpolation='none',cmap=cm,vmin=vmin,vmax=vmax)
 ax.set_yticks(np.arange(len(a)));ax.tick_params(length=0)
 return im
ax=fig.add_axes([.09,.68,.49,.20]);im=heat(ax,x,'viridis',0,1)
ax.set_yticklabels(drugs);ax.set_xticks(range(48),ids,rotation=90,fontsize=7.5)
ax.set_title('A1  PLIF interaction pattern',loc='left',pad=28,fontsize=14)
ax.axvline(20.5,color='white',lw=2)
for pos,lab in [(10,'α1 · F01–F21'),(34,'α2 · F22–F48')]:ax.text(pos,1.02,lab,transform=ax.get_xaxis_transform(),ha='center',fontsize=11)
cb=fig.colorbar(im,cax=fig.add_axes([.20,.635,.26,.008]),orientation='horizontal');cb.set_label('Frozen contact frequency')
ax=fig.add_axes([.61,.68,.12,.20]);im=heat(ax,ranks,'YlGnBu',1,10)
ax.set_yticklabels([]);ax.set_xticks([0,1],['α1','α2']);ax.set_title('A2  Vina score',fontsize=14,pad=28)
for i in range(10):
 for j in range(2):ax.text(j,i,f'{v.iloc[i,j]:.3f}',ha='center',va='center',fontsize=10,color='white' if ranks.iloc[i,j]>6 else '#171717')
cb=fig.colorbar(im,cax=fig.add_axes([.61,.635,.12,.008]),orientation='horizontal',ticks=[1,5,10]);cb.set_label('Within-receptor rank',fontsize=10)
fig.text(.67,.602,'10 = most favorable\nCell text: raw median score (kcal/mol)',ha='center',fontsize=9)
ax=fig.add_axes([.77,.68,.20,.20]);im=heat(ax,y,'YlOrRd',0,float(np.nanmax(y)))
ax.set_yticklabels([]);ax.set_xticks(range(4),terms,rotation=30,ha='right',fontsize=9)
ax.set_title('A3  Clinical C1 fingerprint',fontsize=14,pad=28)
for i,j in zip(*np.where(y.isna())):ax.text(j,i,'NA',ha='center',va='center',fontsize=9)
cb=fig.colorbar(im,cax=fig.add_axes([.79,.61,.16,.008]),orientation='horizontal');cb.set_label('Frozen logROR · gray = NA',fontsize=10)
fig.text(.07,.565,'B  Existing drug-level correspondence',fontsize=18,weight='bold')
ax=fig.add_axes([.09,.115,.35,.42]);heat(ax,pr,'RdBu_r',-1,1)
ax.set_title('B1  PLIF × C1',loc='left',pad=12,fontsize=15)
ax.set_yticklabels([f'{i} | '+('α1' if k<21 else 'α2') for k,i in enumerate(ids)],fontsize=8)
ax.set_xticks(range(4),terms,rotation=20,ha='right',fontsize=10)
ax.axhline(20.5,color='#222',lw=2)
for i in range(48):
 for j in range(4):ax.text(j,i,str(int(pn.iloc[i,j])),ha='center',va='center',fontsize=7.5,color='white' if abs(pr.iloc[i,j])>.58 else '#222')
ax=fig.add_axes([.58,.425,.35,.10]);heat(ax,vr,'RdBu_r',-1,1)
ax.set_title('B2  Vina favorability × C1',loc='left',pad=12,fontsize=15)
ax.set_yticklabels(['α1','α2'],fontsize=12);ax.set_xticks(range(4),terms,rotation=20,ha='right',fontsize=10)
for i in range(2):
 for j in range(4):ax.text(j,i,str(int(vn.iloc[i,j])),ha='center',va='center',fontsize=13,color='white' if abs(vr.iloc[i,j])>.58 else '#222')
cb=fig.colorbar(ScalarMappable(norm=Normalize(-1,1),cmap='RdBu_r'),cax=fig.add_axes([.60,.35,.30,.012]),orientation='horizontal',ticks=[-1,-.5,0,.5,1]);cb.set_label('Spearman rho · shared scale for B1 and B2')
fig.text(.58,.30,'B1 / B2 cell text = n_effective (drugs)\nGray = non-estimable / NA\nAll rho values are read from existing tables.',fontsize=12,linespacing=1.7,va='top')
fig.text(.58,.225,'Views remain independent\n\nPLIF: frozen interaction-pattern frequencies\nVina: median of five independent seed-best scores\nClinical: four prespecified C1 logROR terms\n\nNo composite score or clinical optimization.',fontsize=12,linespacing=1.6,va='top')
fig.text(.07,.063,'PLIF and Vina score are not combined into a composite predictor.\nVina score represents docking-score favorability, not experimental affinity.\nAll correspondence analyses use drug as the observational unit.\nAssociation/correspondence only; not causation or mediation.',fontsize=11,linespacing=1.5,va='top')
for ext in ['png','svg']:fig.savefig(F/f'multiview_plif_vina_c1_profile.{ext}',dpi=220)
plt.close(fig)
# Verify snapshots remain byte-identical after plotting.
for m in manifest:assert hashlib.sha256((R/m['snapshot']).read_bytes()).hexdigest()==m['sha256']
(R/'QC_REPORT.json').write_text(json.dumps({'drug_order_verified':True,'frozen_feature_order_verified':True,'plif_shape':list(x.shape),'clinical_shape':list(y.shape),'plif_correspondence_shape':list(pr.shape),'vina_correspondence_shape':list(vr.shape),'source_hashes_unchanged':True,'correlations_recomputed':False,'display_transform':'within-receptor average rank of negative median score; 1 least, 10 most favorable','clinical_NA_cells':int(y.isna().sum().sum())},indent=2)+'\n')
print('Created PNG/SVG; source hashes and frozen orders verified.')
