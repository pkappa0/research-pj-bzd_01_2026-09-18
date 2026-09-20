from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
R=Path(__file__).resolve().parents[1];D=['diazepam','triazolam','brotizolam','lormetazepam'];T=['Ataxia','Balance disorder','Coordination abnormal','Gait disturbance']
x=pd.read_csv(R/'data/nishino_primary4_plif.csv').set_index('drug').loc[D]
v=pd.read_csv(R/'data/nishino_primary4_vina_summary.csv').pivot(index='drug',columns='receptor_block',values='median_best_seed_score').loc[D,['alpha1','alpha2']]
m=pd.read_csv(R/'data/nishino2008_primary4_in_vivo_M.csv').set_index('drug').loc[D]
y=pd.read_csv(R/'data/nishino_primary4_clinical_Y.csv').set_index('drug').loc[D]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none'})
def save(fig,name):
 for e in ['png','svg']:fig.savefig(R/'figures'/f'{name}.{e}',dpi=220)
 plt.close(fig)
def heat(ax,a,cmap='viridis',vmin=None,vmax=None):
 cm=plt.get_cmap(cmap).copy();cm.set_bad('#c7c7c7');return ax.imshow(np.ma.masked_invalid(np.asarray(a,float)),aspect='auto',cmap=cm,vmin=vmin,vmax=vmax,interpolation='none')
fig=plt.figure(figsize=(26,8))
fig.suptitle('Nishino 2008 matched model case:\nmulti-view structural, in vivo, and clinical phenotype representation',fontsize=22,y=.98)
fig.text(.5,.845,'Four direct-acting benzodiazepines aligned by drug identity',ha='center',fontsize=16)
ax=fig.add_axes([.075,.35,.425,.40]);im=heat(ax,x,vmin=0,vmax=1);ax.set_yticks(range(4),D);ax.set_xticks(range(48),[f'F{i:02}' for i in range(1,49)],rotation=90,fontsize=7)
ax.set_title('A  Structural PLIF',loc='left',pad=30);ax.axvline(20.5,color='white',lw=2)
for pos,label in [(10,'α1 · F01–F21'),(34,'α2 · F22–F48')]:ax.text(pos,1.025,label,transform=ax.get_xaxis_transform(),ha='center')
fig.colorbar(im,cax=fig.add_axes([.16,.23,.24,.014]),orientation='horizontal',label='Contact frequency')
ax=fig.add_axes([.53,.35,.09,.40]);r=(-v).rank();im=heat(ax,r,'YlGnBu',1,4);ax.set_title('B  Vina docking score',pad=30);ax.set_yticks([]);ax.set_xticks([0,1],['α1','α2'])
for i in range(4):
 for j in range(2):ax.text(j,i,f'{v.iloc[i,j]:.3f}',ha='center',va='center',fontsize=11,color='white' if r.iloc[i,j]>2.5 else 'black')
fig.colorbar(im,cax=fig.add_axes([.53,.23,.09,.014]),orientation='horizontal',ticks=[1,2,3,4]);fig.text(.575,.14,'Within-receptor favorability rank\n4 = most favorable; cell = raw kcal/mol',ha='center',fontsize=10)
ax=fig.add_axes([.66,.35,.12,.40]);im=heat(ax,m[['rotarod_potency']],'PuBu',float(m.rotarod_potency.min()),float(m.rotarod_potency.max()));ax.set_title('C  Rotarod phenotype',pad=30);ax.set_xticks([]);ax.set_yticks([])
for i,d in enumerate(D):
 q=m.loc[d];ax.text(0,i,f'ED50 {q.ed50_mg_kg:.2f} mg/kg\n[{q.ed50_ci_low:.2f}, {q.ed50_ci_high:.2f}]\nPotency {q.rotarod_potency:+.3f}',ha='center',va='center',fontsize=10,color='white' if q.rotarod_potency>-.4 else 'black')
fig.colorbar(im,cax=fig.add_axes([.66,.23,.12,.014]),orientation='horizontal');fig.text(.72,.14,'Color: −log10(ED50 mg/kg)\nBrackets: reported interval*',ha='center',fontsize=10)
ax=fig.add_axes([.82,.35,.165,.40]);im=heat(ax,y,'RdBu_r',-float(np.nanmax(np.abs(y))),float(np.nanmax(np.abs(y))));ax.set_title('D  Clinical C1 fingerprint',pad=30);ax.set_yticks([]);ax.set_xticks(range(4),T,rotation=25,ha='right',fontsize=9)
for i,j in zip(*np.where(y.isna())):ax.text(j,i,'NA',ha='center',va='center')
fig.colorbar(im,cax=fig.add_axes([.85,.19,.12,.014]),orientation='horizontal',label='logROR · gray = unavailable')
fig.text(.075,.035,'Representation / hypothesis-generation only. No causal, mediation, or predictive inference.\n*The paper does not explicitly state interval confidence level or ED50 time aggregation. Vina docking score is not experimental affinity.',fontsize=11)
save(fig,'nishino_primary4_XMY_model_case')
raw=pd.read_csv(R/'data/nishino2008_rotarod_raw.csv')
fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(3,2,height_ratios=[1,1,.85],hspace=.45,wspace=.25)
for d,slot in zip(D,[(0,0),(0,1),(1,0),(1,1)]):
 ax=fig.add_subplot(gs[slot]);sub=raw[raw.drug==d]
 for dose,g in sub.groupby('dose_mg_kg'):ax.plot(g.time_min,g.failure_fraction,'-o',label=f'{dose:g} mg/kg')
 ax.set(title=d,xlabel='Time after oral dose (min)',ylabel='Failure fraction',xticks=[15,30,60,90],ylim=(-.03,1.05));ax.grid(alpha=.15);ax.legend(ncol=2,fontsize=9)
ax=fig.add_subplot(gs[2,:]);ax.errorbar(m.ed50_mg_kg,np.arange(4),xerr=np.array([m.ed50_mg_kg-m.ed50_ci_low,m.ed50_ci_high-m.ed50_mg_kg]),fmt='o',color='#246478',capsize=4)
ax.set(yticks=range(4),yticklabels=D,xlabel='Reported rotarod ED50 (mg/kg); horizontal bars = reported intervals*',xlim=(0,11.5));ax.invert_yaxis()
for i,d in enumerate(D):q=m.loc[d];ax.text(q.ed50_ci_high+.2,i,f'{q.ed50_mg_kg:g} [{q.ed50_ci_low:g}, {q.ed50_ci_high:g}]',va='center',fontsize=10)
fig.suptitle('Nishino 2008 rotarod phenotype · primary four drugs',fontsize=19,y=.99)
fig.subplots_adjust(bottom=.12,top=.94)
fig.text(.08,.025,'Male ICR mice; 15 rpm; fall within 180 s; n=10 per dose/time cell. Doses are not normalized across drugs.\nTriazolam 0.5 mg/kg at 90 min = 9/10 as printed. No correction or outlier exclusion.\n*ED50 interval confidence level and time aggregation are not explicitly specified; reported values are not refitted.',fontsize=10)
save(fig,'nishino_primary4_rotarod_profiles')
a=pd.read_csv(R/'data/nishino_primary4_X_M_correspondence.csv');fig,(ax,bx)=plt.subplots(1,2,figsize=(12,15),gridspec_kw={'width_ratios':[1.4,1]})
q=a.iloc[:48];im=heat(ax,q[['spearman_rho']],'RdBu_r',-1,1);ax.set_yticks(range(48),[f'F{i:02} | '+('α1' if i<=21 else 'α2') for i in range(1,49)],fontsize=9);ax.set_xticks([0],['PLIF ↔ rotarod potency']);ax.axhline(20.5,color='black',lw=1)
for i,row in q.iterrows():ax.text(0,i,'NA' if pd.isna(row.spearman_rho) else f'{row.spearman_rho:+.2f}',ha='center',va='center',fontsize=8,color='white' if abs(row.spearman_rho)>.6 else 'black')
bx.axis('off');text='Vina favorability ↔ rotarod potency\n\n'+'\n'.join(f"{r.structural_feature.split('_')[0]}: rho = {r.spearman_rho:+.2f}; n = {r.n_effective}" for r in a.iloc[48:].itertuples())+'\n\nAll PLIF features in frozen order.\nGray = constant / non-estimable.\n\nDrug is the observational unit.\nOnly four drugs; no p-value testing,\nfeature selection or discovery claims.';bx.text(.02,.9,text,va='top',fontsize=12,linespacing=1.8)
fig.colorbar(im,ax=bx,orientation='horizontal',fraction=.05,pad=.2,label='Spearman rho',ticks=[-1,0,1]);fig.suptitle('Descriptive correspondence in a four-drug matched model case\nStructural X ↔ rotarod motor-impairment potency · n=4',fontsize=17);fig.subplots_adjust(top=.93,bottom=.06,left=.13,wspace=.25)
save(fig,'nishino_primary4_X_M_correspondence')
fig,(ax,bx)=plt.subplots(1,2,figsize=(16,7),gridspec_kw={'width_ratios':[1.2,1]});ax.axis('off')
for yy,label in [(.8,'Structural X\n[PLIF α1 | PLIF α2 | Vina α1 | Vina α2]'),(.49,'in vivo M\n[Nishino rotarod ED50]'),(.18,'Clinical Y\n[Balance / Ataxia-related MedDRA PTs]')]:ax.text(.5,yy,label,ha='center',va='center',fontsize=15,bbox=dict(boxstyle='round,pad=.8',facecolor='#eaf2f5',edgecolor='#507789'),transform=ax.transAxes)
for lo,hi in [(.59,.69),(.28,.38)]:ax.plot([.5,.5],[lo,hi],ls=':',color='#64747b',transform=ax.transAxes)
ax.text(.04,.98,'Aligned by drug identity · no causal arrows',transform=ax.transAxes,fontsize=12)
master=pd.read_csv(R/'data/nishino_primary4_XMY_master.csv').set_index('drug').loc[D]
status=np.array([[2 if master.loc[d,'structural_completeness']=='complete' else 1,2 if np.isfinite(m.loc[d,'ed50_mg_kg']) else 0,2 if y.loc[d].notna().sum()==4 else 1 if y.loc[d].notna().any() else 0] for d in D])
im=bx.imshow(status,cmap=ListedColormap(['#bdbdbd','#efc577','#62a999']),vmin=0,vmax=2,aspect='auto');bx.set_xticks(range(3),['Structural X','in vivo M*','Clinical Y']);bx.set_yticks(range(4),D)
for i in range(4):
 for j in range(3):bx.text(j,i,('Complete\n4/4 PTs' if j==2 and status[i,j]==2 else 'Partial\n3/4 PTs' if j==2 else 'Complete'),ha='center',va='center',fontsize=12)
fig.suptitle('Nishino 2008 · matched three-layer representation',fontsize=21);fig.subplots_adjust(top=.84,bottom=.22,wspace=.30)
fig.text(.06,.055,'*M: reported ED50 and intervals are available; interval confidence level and ED50 time aggregation remain unspecified.\nClinical gaps: Coordination abnormal for triazolam and brotizolam. No imputation.\nRepresentation / hypothesis-generation only. No causal, mediation, or predictive inference.',fontsize=11)
save(fig,'nishino_three_layer_summary')
print('Generated four figure sets.')
