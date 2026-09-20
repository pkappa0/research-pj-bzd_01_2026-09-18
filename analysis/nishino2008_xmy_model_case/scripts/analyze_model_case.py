"""Descriptive four-drug correspondence; no learning or feature selection."""
from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1];D=['diazepam','triazolam','brotizolam','lormetazepam']
terms=['Ataxia','Balance disorder','Coordination abnormal','Gait disturbance'];yc=['logROR_ataxia','logROR_balance_disorder','logROR_coordination_abnormal','logROR_gait_disturbance']
cfg={'observational_unit':'drug','n_primary_drugs':4,'order':D,'minimum_pairs':3,'primary_complete_set':4,'partial_set_policy':'n=3 descriptive only, explicitly labeled; n<3 NA','constant_input':'NA','p_values':'not calculated/reported','feature_order':'frozen F01-F48, no ranking','C1_terms':terms,'note':'User-authorized n=4 model case uses n>=3 descriptive floor; historical 10-drug n>=5 results unchanged. No inferential claims.'}
(R/'config/correspondence_config.json').write_text(json.dumps(cfg,indent=2)+'\n')
freeze=json.loads((R/'config/STRUCTURAL_RESULTS_FROZEN.json').read_text())
for name,h in freeze['files'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==h
x=pd.read_csv(R/'data/nishino_primary4_plif.csv').set_index('drug').loc[D]
s=pd.read_csv(R/'data/nishino_primary4_structural_multiview.csv').set_index('drug').loc[D]
m=pd.read_csv(R/'data/nishino2008_primary4_in_vivo_M.csv').set_index('drug').loc[D]
y=pd.read_csv(R/'data/nishino_primary4_clinical_Y.csv').set_index('drug').loc[D]
master=s.join(m[['ed50_mg_kg','ed50_ci_low','ed50_ci_high','rotarod_potency','interval_level','ed50_time_basis']]).join(y)
master['structural_completeness']=['complete' if x.loc[d].notna().all() and s.loc[d].notna().all() else 'partial' for d in D]
master['M_completeness']='complete_reported_values; ED50_time_basis_and_interval_level_unspecified'
master['Y_n_available']=y.notna().sum(axis=1);master['Y_completeness']=np.where(master.Y_n_available==4,'complete','partial')
master.to_csv(R/'data/nishino_primary4_XMY_master.csv',na_rep='NA')
features={f'F{i+1:02}':x[c] for i,c in enumerate(x)}
features.update({b+'_Vina_favorability':s[b+'_Vina_favorability'] for b in ['alpha1','alpha2']})
def corr(a,b):
 mask=np.isfinite(a)&np.isfinite(b);aa=a[mask];bb=b[mask];n=int(mask.sum());status='complete_four_drug' if n==4 else 'partial_drug_set'
 if n<3:rho=np.nan;status='insufficient_clinical_coverage_n_lt_3'
 elif aa.nunique()<2 or bb.nunique()<2:rho=np.nan;status='constant_input'
 else:rho=float(spearmanr(aa,bb).statistic)
 return dict(spearman_rho=rho,n_effective=n,status=status,included_drugs=';'.join(a.index[mask]),excluded_drugs=';'.join(a.index[~mask]))
rows=[]
for fid,a in features.items():rows.append(dict(structural_feature=fid,feature_label=x.columns[int(fid[1:])-1] if fid.startswith('F') else fid,outcome='rotarod_potency',**corr(a,m.rotarod_potency)))
pd.DataFrame(rows).to_csv(R/'data/nishino_primary4_X_M_correspondence.csv',index=False,na_rep='NA')
rows=[dict(clinical_term=t,**corr(m.rotarod_potency,y[c])) for t,c in zip(terms,yc)]
pd.DataFrame(rows).to_csv(R/'data/nishino_primary4_M_Y_correspondence.csv',index=False,na_rep='NA')
rows=[dict(structural_feature=fid,clinical_term=t,secondary=True,**corr(a,y[c])) for fid,a in features.items() for t,c in zip(terms,yc)]
pd.DataFrame(rows).to_csv(R/'data/nishino_primary4_X_Y_correspondence.csv',index=False,na_rep='NA')
summary={'X_complete_drugs':int((master.structural_completeness=='complete').sum()),'M_reported_ED50_complete_drugs':int(m.ed50_mg_kg.notna().sum()),'Y_available_cells':int(y.notna().sum().sum()),'Y_total_cells':16,'Y_complete_drugs':int(y.notna().all(axis=1).sum()),'X_variable_features':int((x.nunique()>1).sum()),'X_constant_features':int((x.nunique()==1).sum()),'ED50_range':[float(m.ed50_mg_kg.min()),float(m.ed50_mg_kg.max())],'M_dynamic_range_fold':float(m.ed50_mg_kg.max()/m.ed50_mg_kg.min()),'no_ML':True,'no_composite':True,'no_inference':True,'decision':'CASE1 qualified: 3/4 PTs complete across four drugs; targeted CASE2 coverage gap for Coordination abnormal. M summary-method uncertainty remains.'}
(R/'qc/model_case_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2));print(pd.read_csv(R/'data/nishino_primary4_M_Y_correspondence.csv').to_string(index=False))
