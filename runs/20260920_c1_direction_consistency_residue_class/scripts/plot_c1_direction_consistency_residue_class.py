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
SOURCE=ROOT/'runs/20260920_c1_direction_consistency'
DATA=RUN/'data';FIG=RUN/'figures'
DATA.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def annotation(n,pos,neg):
 if n==0:return 'NA (0/4 estimable)'
 if n<4:return f'mixed ({n}/4 estimable)'
 if pos==4:return '4/4 positive'
 if neg==4:return '4/4 negative'
 return 'mixed'
def compact(label,fid):
 receptor,site,residue,kind=label.split('|');assert kind=='any_contact_frequency'
 m=re.fullmatch(r'BZD_(GAMMA2|SITE)_(\d+)',site);assert m
 return f"{fid} | {'α1' if receptor=='alpha1' else 'α2'} | {'γ2' if m[1]=='GAMMA2' else 'site'} {residue}{int(m[2])}"
manifest=[]
for name,src in [('c1_direction_consistency_summary.csv',SOURCE/'data/c1_direction_consistency_summary.csv'),('feature_display_key.csv',SOURCE/'data/feature_display_key.csv'),('ANALYSIS_CONTRACT_source.md',ROOT/'ANALYSIS_CONTRACT.md')]:
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
classes=[('Aromatic',['PHE','TYR','TRP']),('Hydrophobic / non-aromatic',['ALA','VAL','ILE','LEU','MET','PRO','GLY']),('Hydrophilic / polar uncharged',['SER','THR','ASN','GLN','CYS']),('Acidic / ionizable',['ASP','GLU']),('Basic / ionizable',['LYS','ARG','HIS'])]
classmap={res:label for label,members in classes for res in members}
s['residue']=s.feature_label.str.split('|').str[2]
assert s.residue.isin(classmap).all()
s['residue_class']=s.residue.map(classmap)
s['class_order']=s.residue_class.map({label:i for i,(label,_) in enumerate(classes)})
s=s.sort_values(['class_order','median_rho_across_c1_terms'],ascending=[True,False],kind='stable',na_position='last').reset_index(drop=True)
headers=[];positions=[];cursor=0
for label,members in classes:
 count=int((s.residue_class==label).sum())
 if count:
  headers.append((cursor,label,count));positions.extend(range(cursor+1,cursor+1+count));cursor+=count+2
s['plot_y']=positions
mapping=s[['feature_id','feature_label','residue','residue_class']].copy()
mapping['aromatic_sidechain_flag']=mapping.residue.isin(['PHE','TYR','TRP','HIS'])
mapping['polar_sidechain_note']=mapping.residue.map({'TYR':'phenolic hydroxyl; aromatic and polar','TRP':'indole nitrogen; aromatic','HIS':'aromatic imidazole; placed in basic/ionizable display group'}).fillna('')
mapping.to_csv(DATA/'residue_class_mapping.csv',index=False)
(RUN/'residue_class_rules.json').write_text(json.dumps({'display_groups':dict(classes),'priority_notes':['PHE/TYR/TRP in aromatic group despite overlapping polarity','HIS in basic/ionizable group; aromatic flag retained','Ionization labels do not infer actual protonation in receptor','No hydropathy score calculated'],'within_group_order':'median rho descending, stable frozen-order ties'},indent=2)+'\n')
labels=[compact(f,fid) for fid,f in zip(s.feature_id,s.feature_label)]
notes=[annotation(int(r.number_estimable_c1_terms),int(r.number_positive_rho),int(r.number_negative_rho)) for r in s.itertuples()]
# Additional metadata only; values in the input summary are never recalculated.
audit=s[['feature_id','median_rho_across_c1_terms','min_rho','max_rho']].copy();audit['compact_label']=labels;audit['right_annotation']=notes;audit['residue_class']=s.residue_class.values
audit.to_csv(DATA/'figure_display_annotations.csv',index=False,na_rep='NA')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False})
fig=plt.figure(figsize=(13.5,23),facecolor='white');ax=fig.add_axes([.30,.13,.47,.755]);ann=fig.add_axes([.80,.13,.19,.755],sharey=ax);ann.axis('off')
fig.text(.5,.975,'C1 Balance / Ataxia: correspondence grouped by residue class',ha='center',fontsize=18,fontweight='bold')
fig.text(.5,.951,'Median and range of drug-level Spearman rho across four prespecified C1 terms',ha='center',fontsize=12.5)
fig.text(.30,.925,'●  Median     ━  Min–max across clinical terms',fontsize=12,color='#216778')
fig.text(.80,.911,'Direction across C1 terms',fontsize=11,fontweight='bold')
fig.text(.30,.897,'● α1',fontsize=11,fontweight='bold',color='#236f89')
fig.text(.40,.897,'● α2',fontsize=11,fontweight='bold',color='#ad6329')
colors={'alpha1':'#236f89','alpha2':'#ad6329'}
ax.axvline(0,color='#65717d',ls='--',lw=1.1,zorder=1)
for i in range(48):
 yy=s.iloc[i].plot_y
 if i%2==0:
  ax.axhspan(yy-.5,yy+.5,color='#f5f7f9',zorder=0);ann.axhspan(yy-.5,yy+.5,color='#f5f7f9',zorder=0)
 r=s.iloc[i]
 color=colors[r.feature_label.split('|')[0]]
 if r.number_estimable_c1_terms>0:
  ax.hlines(yy,r.min_rho,r.max_rho,color=color,lw=2.2,zorder=2)
  ax.scatter(r.median_rho_across_c1_terms,yy,s=35,color=color,edgecolors='white',linewidth=.5,zorder=3)
 else:ax.text(0,yy,'NA',ha='center',va='center',color='#777777',fontsize=10,bbox={'facecolor':'white','edgecolor':'none','pad':1})
 ann.text(.015,yy,notes[i],ha='left',va='center',fontsize=10.5,color='#37414c')
for header,label,count in headers:
 for a in [ax,ann]:a.axhspan(header-.38,header+.38,color='#e8edf1',zorder=0)
 ax.text(-2.04,header,f'{label}  ({count} features)',va='center',fontsize=12,fontweight='bold',color='#354453',clip_on=False)
for a in [ax,ann]:a.set_ylim(cursor-1.3,-.8)
ax.set_xlim(-1,1);ax.set_xticks([ -1,-.75,-.5,-.25,0,.25,.5,.75,1]);ax.set_xlabel('Spearman rho',fontsize=13,labelpad=12)
ax.set_yticks(positions,labels,fontsize=10.5);ax.tick_params(axis='y',length=0,pad=11);ax.tick_params(axis='x',labelsize=11)
ax.grid(axis='x',color='#e3e8ec',lw=.5,zorder=0)
fig.text(.30,.074,'Grouped by residue class; median rho descending within groups. Descriptive only; no feature selection.\n4/4 = all four terms share the sign; mixed = otherwise. Incomplete coverage is annotated.',fontsize=10.5,linespacing=1.55)
fig.text(.07,.024,'Observational unit = drug.\nIntervals show variation across clinical terms, not confidence intervals.\nAssociation/correspondence only; no causal or mechanistic inference.',fontsize=12,linespacing=1.5)
for ext in ['png','svg']:fig.savefig(FIG/f'c1_plif_direction_consistency_residue_class.{ext}',dpi=240,facecolor='white')
plt.close(fig)
qc={'input_summary_sha256':sha(DATA/'c1_direction_consistency_summary.csv'),'feature_key_sha256':sha(DATA/'feature_display_key.csv'),'feature_count':len(s),'feature_order':s.feature_id.tolist(),'alpha1_feature_count':21,'alpha2_feature_count':27,'estimable_rows':int(present.sum()),'undefined_rows':int((~present).sum()),'axis_limits':[-1,1],'correlations_recomputed':False,'summary_statistics_recomputed':False,'feature_selection':False,'display_sort':'residue class then median descending; stable frozen-order ties; NA last','class_counts':s.residue_class.value_counts().to_dict(),'interval_meaning':'min/max rho across C1 terms, not confidence intervals','annotations':dict(zip(s.feature_id,notes))}
(RUN/'QC_REPORT.json').write_text(json.dumps(qc,indent=2,ensure_ascii=False)+'\n')
print('Created PNG/SVG. All 48 rows retained; existing medians/ranges used verbatim.')
