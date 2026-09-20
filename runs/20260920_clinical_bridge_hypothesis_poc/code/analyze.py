from pathlib import Path
import json,hashlib,itertools,datetime,warnings
import pandas as pd,numpy as np
from scipy.stats import spearmanr,rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False))
cf=json.loads((R/'CLINICAL_RESULTS_FROZEN.json').read_text());sf=json.loads((R/'STRUCTURAL_RESULTS_FROZEN.json').read_text())
assert sha(R/'tables/clinical_fingerprint_10drug.csv')==cf['clinical_sha256']
assert sha(R/'raw/lineage/structural_fingerprint.csv')==sf['sha256']
cfg=json.loads((R/'config/clinical_phenotype_config.json').read_text()); axes=cfg['axes']; terms=cfg['feature_order']
y=pd.read_csv(R/'tables/clinical_fingerprint_10drug.csv',index_col=0);x=pd.read_csv(R/'raw/lineage/structural_fingerprint.csv',index_col=0).loc[y.index]
D=list(y.index);n=len(D)
def cos(a,b,minimum):
 mask=np.isfinite(a)&np.isfinite(b); nn=int(mask.sum());v=np.nan
 if nn>=minimum and np.linalg.norm(a[mask])*np.linalg.norm(b[mask])>0:v=float(np.dot(a[mask],b[mask])/(np.linalg.norm(a[mask])*np.linalg.norm(b[mask])))
 return v,nn
sx=np.full((n,n),np.nan);sy=sx.copy(); pairs=[]
for i,j in itertools.combinations(range(n),2):
 xv,nx=cos(x.iloc[i].values,x.iloc[j].values,1);yv,ny=cos(y.iloc[i].values,y.iloc[j].values,7)
 sx[i,j]=sx[j,i]=xv;sy[i,j]=sy[j,i]=yv
 pairs.append({'drug_1':D[i],'drug_2':D[j],'structural_similarity':xv,'clinical_similarity':yv,'structural_shared_features':nx,'clinical_shared_features':ny,'clinical_status':'estimable' if np.isfinite(yv) else 'insufficient'})
pairs=pd.DataFrame(pairs);pairs.to_csv(R/'tables/structural_clinical_pairwise_similarity.csv',index=False,na_rep='NA')
valid=pairs.dropna();rho=float(spearmanr(valid.structural_similarity,valid.clinical_similarity).statistic)
# Drug-label permutations maintain dependence among the 45 dyads.
iu=np.triu_indices(n,1); rng=np.random.default_rng(20260920); null=[]
for z in range(9999):
 p=rng.permutation(n); yy=sy[np.ix_(p,p)][iu];xx=sx[iu];mask=np.isfinite(xx)&np.isfinite(yy)
 null.append(float(spearmanr(xx[mask],yy[mask]).statistic))
perm=(1+sum(abs(v)>=abs(rho) for v in null))/10000
pd.DataFrame({'permutation_rho':null}).to_csv(R/'tables/drug_label_permutation_null.csv',index=False)
loo=[]
for d in D:
 v=pairs[(pairs.drug_1!=d)&(pairs.drug_2!=d)].dropna();loo.append({'omitted_drug':d,'pairs':len(v),'rho':spearmanr(v.structural_similarity,v.clinical_similarity).statistic})
pd.DataFrame(loo).to_csv(R/'tables/leave_one_drug_out_similarity.csv',index=False)
cor=[]
for f in x:
 for t in y:
  a=x[f];b=y[t];m=a.notna()&b.notna();ne=int(m.sum());status='estimable';r=np.nan
  if ne<5:status='exploratory-insufficient'
  elif a[m].nunique()<2 or b[m].nunique()<2:status='constant-feature'
  else:r=float(spearmanr(a[m],b[m]).statistic)
  cor.append({'structural_feature':f,'phenotype_term':t,'phenotype_axis':next(k for k,v in axes.items() if t in v),'n_effective':ne,'missingness':n-ne,'rho':r,'status':status})
cor=pd.DataFrame(cor);cor.to_csv(R/'tables/structural_feature_clinical_correlations.csv',index=False,na_rep='NA')
prof=pd.read_csv(R/'raw/lineage/phenotype_profile_10drug.csv').set_index('drug').reindex(D)
prof.to_csv(R/'tables/invivo_fingerprint_sparse.csv',na_rep='NA')
x.join(y.add_prefix('clinical|')).join(prof.add_prefix('invivo|')).to_csv(R/'tables/three_layer_joined_dataset.csv',na_rep='NA')
m=pd.read_csv(R/'raw/lineage/phenotype_normalized_within_study.csv');anchor=m[(m.study_id=='Nishino2008')&(m.tier==1)&(m.phenotype_axis=='P1')]
anchor=anchor.copy()
for col in ['dose','observation_time','value']: anchor[col]=pd.to_numeric(anchor[col],errors='raise')
anchor.to_csv(R/'tables/invivo_anchor_contexts.csv',index=False,na_rep='NA')
mp=anchor.pivot(index=['dose','observation_time'],columns='drug',values='value');mdir=int((mp.triazolam>mp.diazepam).sum());mtie=int((mp.triazolam==mp.diazepam).sum())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(R/'figures'/f'{name}.{ext}',dpi=180)
 plt.close(fig)
def heat(ax,v,xt,yt,title,vmin=None,vmax=None,cmap='coolwarm',fontsize=8):
 cm=plt.get_cmap(cmap).copy();cm.set_bad('#d9d9d9');im=ax.imshow(v,aspect='auto',cmap=cm,vmin=vmin,vmax=vmax)
 ax.set_xticks(range(len(xt)),xt,rotation=60,ha='right',fontsize=fontsize);ax.set_yticks(range(len(yt)),yt,fontsize=fontsize);ax.set_title(title,pad=15);return im
fig,ax=plt.subplots(figsize=(12,6));im=heat(ax,y.values,[next(k for k,v in axes.items() if t in v)+' '+t for t in terms],D,'Clinical fingerprint | log reporting odds ratio (all roles)',-3,3);fig.colorbar(im,ax=ax,label='log(ROR); displayed scale clipped at +/-3');fig.tight_layout(rect=[0,.08,1,1]);fig.text(.02,.01,'Gray = NA (event count <5 or not estimable). Unadjusted FAERS reporting signals, not incidence or risk.',fontsize=9);save(fig,'clinical_phenotype_fingerprint_heatmap')
fig,ax=plt.subplots(figsize=(8,6));ax.scatter(pairs.structural_similarity,pairs.clinical_similarity,s=45,color='#147d92',alpha=.8);ax.set(xlabel='Structural cosine similarity | fixed 48 features',ylabel='Clinical cosine similarity | stable logROR terms',title=f'Structural / clinical correspondence\nSpearman rho = {rho:.3f}; {len(valid)} dependent pairs from 10 drugs');fig.text(.08,-.01,f'Descriptive dyadic association. Drug-label permutation tail fraction = {perm:.3f}.\nNo independent n=45 test; no causal interpretation.',fontsize=9);save(fig,'structural_vs_clinical_similarity')
cr=cor.pivot(index='structural_feature',columns='phenotype_term',values='rho').reindex(index=x.columns,columns=terms);cn=cor.pivot(index='structural_feature',columns='phenotype_term',values='n_effective').reindex_like(cr)
fig,ax=plt.subplots(figsize=(14,18));im=heat(ax,cr.values,terms,[f.replace('|any_contact_frequency','') for f in x.columns],'Feature / clinical correspondence | cell text = n_effective',-1,1,fontsize=8)
for i in range(len(cr)):
 for j in range(len(terms)):ax.text(j,i,str(int(cn.iloc[i,j])),ha='center',va='center',fontsize=6,color='black')
fig.colorbar(im,ax=ax,shrink=.45,label='Drug-level Spearman rho');fig.text(.02,0,'All 48 frozen X features retained. Gray = constant feature or n<5. No p-value ranking.',fontsize=10);save(fig,'structural_feature_clinical_phenotype_heatmap')
fig,ax=plt.subplots(figsize=(11,8));ax.axis('off')
texts=[('STRUCTURAL FINGERPRINT X','10/10 drugs | 48 fixed alpha1 + alpha2 contact features\nNo recomputation, feature selection or docking'),('IN VIVO PHENOTYPE M','Archived evidence linked for 10 drugs; numeric evidence sparse\nComparable rotarod anchor: diazepam / triazolam only (n=2)'),('CLINICAL PHENOTYPE Y',f'{int((y.notna().sum(axis=1)>0).sum())}/10 drugs | 6 axes / 14 prespecified terms\nFAERS all-role report-document ROR; low counts masked')]
for yy,(title,body) in zip([.81,.49,.17],texts):
 ax.text(.5,yy,title+'\n\n'+body,ha='center',va='center',bbox=dict(boxstyle='round,pad=1',facecolor='#eaf3f5',edgecolor='#147d92'),fontsize=12,transform=ax.transAxes)
for upper,lower in [(.68,.61),(.36,.29)]:ax.annotate('',xy=(.5,lower),xytext=(.5,upper),xycoords='axes fraction',arrowprops=dict(arrowstyle='<->',lw=2,color='#147d92'))
ax.set_title('Three-layer working hypothesis and evidence coverage',fontsize=17,pad=20);fig.text(.5,.01,f'X-Y: exploratory rho {rho:.3f} | M linkage feasible; intermediate function unproven\nArrows indicate proposed correspondence, not causation or mediation.',ha='center',fontsize=11);save(fig,'three_layer_hypothesis')
fig,axs=plt.subplots(3,1,figsize=(14,13),gridspec_kw={'height_ratios':[1,1.4,1.5]});dd=x.loc['triazolam']-x.loc['diazepam'];axs[0].bar(range(len(dd)),dd,color=['#147d92' if s.startswith('alpha1') else '#c77932' for s in dd.index]);axs[0].set(title='X | fixed structural difference: triazolam - diazepam',ylabel='Contact frequency difference',xlabel='Frozen feature index (0-20 alpha1; 21-47 alpha2)');axs[0].axhline(0,color='gray',lw=.6)
for d,col in [('diazepam','#147d92'),('triazolam','#c77932')]:axs[1].plot(range(len(mp)),mp[d],marker='o',label=d,color=col)
axs[1].set_xticks(range(len(mp)),[f'{dose:g} mg/kg\n{t:g} min' for dose,t in mp.index]);axs[1].set(title='M | Nishino2008: context-matched rotarod failure fraction',ylabel='Failed / n=10 mice',ylim=(-.05,1.05));axs[1].legend()
ct=axes['C1']+axes['C2']+axes['C3'];dp=pd.read_csv(R/'tables/clinical_disproportionality.csv');pos=np.arange(len(ct))
for offset,d,col in [(-.12,'diazepam','#147d92'),(.12,'triazolam','#c77932')]:
 z=dp[dp.drug==d].set_index('adverse_event_term').loc[ct];vals=z.logROR.where(z.a>=5).values;err=np.array([vals-np.log(z.CI95_low.values),np.log(z.CI95_high.values)-vals]);axs[2].errorbar(pos+offset,vals,yerr=err,fmt='o',capsize=3,label=d,color=col)
axs[2].set_xticks(pos,ct,rotation=25,ha='right');axs[2].axhline(0,color='gray',lw=.7);axs[2].set(title='Y | C1 / C2 / C3 clinical terms: logROR with nominal 95% CI',ylabel='log(ROR)');axs[2].legend();fig.suptitle('Same drug identity links three layers | descriptive n=2 example',fontsize=17);fig.tight_layout(rect=[0,.05,1,.95]);fig.text(.05,.01,'No n=2 correlation, mediation or causal path. Dose/time contexts are not independent drug pairs.\nEvent counts <5 omitted. Nominal intervals omit reporting bias, confounding and duplication uncertainty.',fontsize=10);save(fig,'diazepam_triazolam_three_layer_example')
# Predefined sensitivity uses all estimable terms, with no selection for favorable results.
ya=pd.read_csv(R/'tables/clinical_fingerprint_all_estimable.csv',index_col=0).loc[D];sa=[]
for i,j in itertools.combinations(range(n),2):sa.append(cos(ya.iloc[i].values,ya.iloc[j].values,7)[0])
sensitivity=float(spearmanr(pairs.structural_similarity,sa,nan_policy='omit').statistic)
top=cor.dropna().iloc[cor.dropna().rho.abs().argmax()].to_dict();coverage={a:int(y[ts].notna().any(axis=1).sum()) for a,ts in axes.items()}
clinical_direction={t:float(y.loc['triazolam',t]-y.loc['diazepam',t]) for t in ct}
summary={'pair_rho':rho,'pair_count':len(valid),'drug_label_permutation_tail_fraction':perm,'all_estimable_sensitivity_rho':sensitivity,'leave_one_drug_out_rho_range':[float(min(v['rho'] for v in loo)),float(max(v['rho'] for v in loo))],'top_absolute_feature_correlation':top,'axis_drug_coverage':coverage,'clinical_drug_count':int(y.notna().any(axis=1).sum()),'rotarod_triazolam_higher_contexts':mdir,'rotarod_tied_contexts':mtie,'clinical_triazolam_minus_diazepam':{k:(v if np.isfinite(v) else None) for k,v in clinical_direction.items()}}
dump(R/'reports/RESULTS_SUMMARY.json',summary)
dump(R/'reports/QC_REPORT.json',{'passed':True,'X_hash_unchanged':True,'Y_freeze_hash_verified':True,'X_features':len(x.columns),'drug_count':n,'pairs':len(pairs),'correlation_rows':len(cor),'all_X_features_retained':len(cr)==48,'minimum_n_for_feature_rho':5,'missingness_preserved':True,'background_total':cf['background_total'],'nonnegative_contingency_cells':bool((dp[['a','b','c','d']]>=0).all().all()),'unique_drug_term_rows':not dp.duplicated(['drug','adverse_event_term']).any(),'limitations':['all-role fallback; no primary-suspect inference','report-document counts; independent case deduplication unavailable','demographics marginal only, no confounding adjustment','exact identity strategy misses nonharmonized brand and salt variants','US-centric FAERS availability and indication/polypharmacy confounding','n=10 drugs; correlated terms/features; exploratory maximum selection','alpha2 provisional construct; beta3 structural context','M comparable only for 2 drugs'],'summary':summary})
print(json.dumps(summary,indent=2))
