#!/usr/bin/env python3
"""Plot archived C1 summaries without computing any correlation or summary statistic."""
from pathlib import Path
import hashlib,json,shutil,re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
RUN=Path(__file__).resolve().parents[1];ROOT=RUN.parents[1]
SOURCE=ROOT/'runs/20260920_wobbling_related_phenotypes'
DATA=RUN/'data';FIG=RUN/'figures'
DATA.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def annotation(n,pos,neg):
 if n<4:return f'{n}/{n} estimable (of 4)'
 if pos==4:return '4/4 +'
 if neg==4:return '4/4 -'
 if pos==3:return '3/4 +'
 if neg==3:return '3/4 -'
 return 'mixed'
def compact(label,fid):
 receptor,site,residue,kind=label.split('|');assert kind=='any_contact_frequency'
 m=re.fullmatch(r'BZD_(GAMMA2|SITE)_(\d+)',site);assert m
 return f"{fid} | {'α1' if receptor=='alpha1' else 'α2'} | {'γ2' if m[1]=='GAMMA2' else 'site'} {residue}{int(m[2])}"
manifest=[]
for name,src in [('c1_direction_consistency_summary.csv',SOURCE/'tables/c1_direction_consistency_summary.csv'),('feature_display_key.csv',SOURCE/'tables/feature_display_key.csv'),('ANALYSIS_CONTRACT_source.md',ROOT/'ANALYSIS_CONTRACT.md')]:
 dest=DATA/name
 if src.exists():
  if dest.exists():assert sha(dest)==sha(src)
  else:shutil.copy2(src,dest)
 elif not dest.exists():raise FileNotFoundError(src)
 manifest.append({'file':name,'source':str(src.relative_to(ROOT)),'sha256':sha(dest)})
if (RUN/'input_manifest.json').exists():
 for v in json.loads((RUN/'input_manifest.json').read_text()):assert sha(DATA/v['file'])==v['sha256']
(RUN/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
s=pd.read_csv(DATA/'c1_direction_consistency_summary.csv');k=pd.read_csv(DATA/'feature_display_key.csv')
ids=[f'F{i:02}' for i in range(1,49)]
assert s.feature_id.tolist()==ids and k.feature_id.tolist()==ids
assert s.feature_label.tolist()==k.feature_label.tolist()
assert s.feature_label.str.startswith('alpha1|').tolist()==[True]*21+[False]*27
assert s.feature_label.iloc[21:].str.startswith('alpha2|').all()
assert (s.number_estimable_c1_terms==s.number_positive_rho+s.number_negative_rho+s.number_zero_rho).all()
assert s.number_estimable_c1_terms.between(0,4).all()
cols=['median_rho_across_c1_terms','min_rho','max_rho'];present=s.number_estimable_c1_terms>0
assert s.loc[present,cols].notna().all().all() and s.loc[~present,cols].isna().all().all()
assert (s.loc[present,'min_rho']<=s.loc[present,'median_rho_across_c1_terms']).all()
assert (s.loc[present,'median_rho_across_c1_terms']<=s.loc[present,'max_rho']).all()
assert (s.loc[present,cols].abs()<=1).all().all()
labels=[compact(f,fid) for fid,f in zip(s.feature_id,s.feature_label)]
notes=[annotation(int(r.number_estimable_c1_terms),int(r.number_positive_rho),int(r.number_negative_rho)) for r in s.itertuples()]
# Additional metadata only; values in the input summary are never recalculated.
audit=s[['feature_id','median_rho_across_c1_terms','min_rho','max_rho']].copy();audit['compact_label']=labels;audit['right_annotation']=notes
audit.to_csv(DATA/'figure_display_annotations.csv',index=False,na_rep='NA')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False})
fig=plt.figure(figsize=(13.5,21),facecolor='white');ax=fig.add_axes([.30,.13,.47,.755]);ann=fig.add_axes([.80,.13,.19,.755],sharey=ax);ann.axis('off')
fig.text(.5,.975,'C1 Balance / Ataxia: consistency of PLIF–clinical correspondence',ha='center',fontsize=18,fontweight='bold')
fig.text(.5,.951,'Median and range of drug-level Spearman rho across four prespecified C1 terms',ha='center',fontsize=12.5)
fig.text(.30,.925,'●  Median     ━  Min–max across clinical terms',fontsize=12,color='#216778')
fig.text(.80,.911,'Direction across C1 terms',fontsize=11,fontweight='bold')
fig.text(.30,.897,'α1 block: F01–F21',fontsize=11,fontweight='bold',color='#414c57')
ax.axvline(0,color='#65717d',ls='--',lw=1.1,zorder=1)
for i in range(48):
 if i%2==0:
  ax.axhspan(i-.5,i+.5,color='#f5f7f9',zorder=0);ann.axhspan(i-.5,i+.5,color='#f5f7f9',zorder=0)
 r=s.iloc[i]
 if r.number_estimable_c1_terms>0:
  ax.hlines(i,r.min_rho,r.max_rho,color='#2f7d8c',lw=2.2,zorder=2)
  ax.scatter(r.median_rho_across_c1_terms,i,s=35,color='#174f61',edgecolors='white',linewidth=.5,zorder=3)
 else:ax.text(0,i,'NA',ha='center',va='center',color='#777777',fontsize=10,bbox={'facecolor':'white','edgecolor':'none','pad':1})
 ann.text(.015,i,notes[i],ha='left',va='center',fontsize=10.5,color='#37414c')
for a in [ax,ann]:a.axhline(20.5,color='#7d8993',lw=1.6,zorder=4);a.set_ylim(47.6,-.8)
ax.text(-1.78,21,'α2',va='center',fontsize=11,fontweight='bold',color='#414c57',clip_on=False)
ax.set_xlim(-1,1);ax.set_xticks([ -1,-.75,-.5,-.25,0,.25,.5,.75,1]);ax.set_xlabel('Spearman rho',fontsize=13,labelpad=12)
ax.set_yticks(range(48),labels,fontsize=10.5);ax.tick_params(axis='y',length=0,pad=11);ax.tick_params(axis='x',labelsize=11)
ax.grid(axis='x',color='#e3e8ec',lw=.5,zorder=0)
fig.text(.30,.074,'+ / − = positive / negative rho; mixed = neither direction has 3 or 4 terms.\nFor <4 estimable terms, the available count is shown explicitly; missing values remain NA.',fontsize=10.5,linespacing=1.55)
fig.text(.07,.024,'Observational unit = drug.\nIntervals show variation across clinical terms, not confidence intervals.\nAssociation/correspondence only; no causal or mechanistic inference.',fontsize=12,linespacing=1.5)
for ext in ['png','svg']:fig.savefig(FIG/f'c1_plif_direction_consistency.{ext}',dpi=240,facecolor='white')
plt.close(fig)
qc={'input_summary_sha256':sha(DATA/'c1_direction_consistency_summary.csv'),'feature_key_sha256':sha(DATA/'feature_display_key.csv'),'feature_count':len(s),'feature_order':ids,'alpha1_feature_count':21,'alpha2_feature_count':27,'estimable_rows':int(present.sum()),'undefined_rows':int((~present).sum()),'axis_limits':[-1,1],'correlations_recomputed':False,'summary_statistics_recomputed':False,'feature_ranking_or_selection':False,'interval_meaning':'min/max rho across C1 terms, not confidence intervals','annotations':dict(zip(ids,notes))}
(RUN/'QC_REPORT.json').write_text(json.dumps(qc,indent=2,ensure_ascii=False)+'\n')
print('Created PNG/SVG. All 48 rows retained; existing medians/ranges used verbatim.')
