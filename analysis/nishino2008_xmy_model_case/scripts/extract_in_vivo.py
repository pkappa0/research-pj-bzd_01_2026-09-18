"""Audited transcription of Nishino 2008 Table 1 and reported Results ED50."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
R=Path(__file__).resolve().parents[1]
# Four entries per dose: 15,30,60,90 min. Includes published unusual 9/10.
raw={
'diazepam':{2:[1,3,4,3],5:[2,4,5,4],10:[4,7,9,8],20:[9,10,10,9]},
'triazolam':{.5:[2,3,3,9],1:[3,4,4,2],2:[5,8,7,5],5:[8,9,9,4]},
'rilmazafone':{5:[1,2,2,2],10:[3,5,5,2],20:[6,7,8,4],30:[8,9,9,6]},
'brotizolam':{2:[0,2,2,1],5:[2,4,4,2],10:[5,6,7,5],20:[6,6,9,5]},
'lormetazepam':{2:[1,4,4,0],5:[2,5,6,2],10:[7,9,9,5],20:[9,10,10,9]}}
ref='https://doi.org/10.1254/jphs.08107FP';rows=[]
for drug,doses in raw.items():
 for dose,counts in doses.items():
  for t,n in zip([15,30,60,90],counts):
   rows.append(dict(study_id='Nishino2008',drug=drug,dose_mg_kg=dose,time_min=t,n_positive=n,n_total=10,failure_fraction=n/10,administration_route='oral',species='Mus musculus',strain='ICR',sex='male',rotarod_rpm=15,trial_duration_sec=180,source_table='Table 1, p353',source_reference=ref,vehicle='0.5% carboxymethylcellulose sodium',administration_volume_ml_kg=10,positive_event='fall within 180 seconds',source_anomaly='published 9/10 without significance marker; retained verbatim' if drug=='triazolam' and dose==.5 and t==90 else ''))
pd.DataFrame(rows).to_csv(R/'data/nishino2008_rotarod_raw.csv',index=False)
pd.DataFrame(rows).pivot(index=['drug','dose_mg_kg'],columns='time_min',values='failure_fraction').to_csv(R/'data/nishino2008_rotarod_dose_time_matrix.csv')
pd.DataFrame([dict(study_id='Nishino2008',drug='vehicle_control',dose_mg_kg=np.nan,time_min=t,n_positive=0,n_total=10,failure_fraction=0,source_table='Table 1, p353',source_reference=ref) for t in [15,30,60,90]]).to_csv(R/'data/nishino2008_rotarod_control.csv',index=False,na_rep='NA')
ed=[]
for d,v,lo,hi in [('diazepam',3.11,2.53,3.81),('triazolam',1.25,1,1.49),('rilmazafone',9.55,8.20,12.9),('brotizolam',5.76,5.45,9.12),('lormetazepam',3.39,2.75,4.27)]:
 ed.append(dict(drug=d,ed50_mg_kg=v,ed50_ci_low=lo,ed50_ci_high=hi,ed50_method='reported probit; not refitted',primary_model_case=d!='rilmazafone',exclusion_reason='prodrug; oral activity reflects active metabolites' if d=='rilmazafone' else '',interval_level='not explicitly specified in article',ed50_time_basis='not explicitly specified; Results discuss 60 min immediately before ED50; not independently assumed',source_reference=ref,source_location='Results pp350-351'))
e=pd.DataFrame(ed);e.to_csv(R/'data/nishino2008_rotarod_ed50.csv',index=False)
m=e[e.primary_model_case].set_index('drug').loc[['diazepam','triazolam','brotizolam','lormetazepam']].copy();m['rotarod_potency']=-np.log10(m.ed50_mg_kg);m.to_csv(R/'data/nishino2008_primary4_in_vivo_M.csv')
protocol={'species':'Mus musculus','strain':'ICR','sex':'male','body_weight_g':'23–28','supplier':'Japan SLC, Shizuoka','route':'oral','vehicle':'0.5% carboxymethylcellulose sodium','administration_volume_ml_kg':10,'rotation_rpm':15,'rod_diameter_cm':3,'compartments':5,'trial_duration_sec':180,'selection':'Animals remaining on the rod in two successive trials 24 hours before the experiment were selected','observation_minutes':[15,30,60,90],'positive_event':'fall within 180 seconds','n_tested_per_dose_time_cell':10,'independent_animal_total':'not explicitly reported; repeated times must not be counted as independent animals','age':'5 weeks stated for plus-maze; not separately specified for rotarod','ed50_method':'probit','ed50_confidence_level':'not explicitly stated','ed50_time_aggregation':'not explicitly stated','source_reference':ref,'source_pages':'349 animals; 350 rotarod/drugs/statistics; 350–351 ED50; 353 Table1'}
(R/'data/nishino2008_protocol.json').write_text(json.dumps(protocol,indent=2,ensure_ascii=False)+'\n')
p=Path('/Users/k-atsumi/Documents/Codex/2026-09-18/github-readme-final-report-csv-git/work/phenotype_poc/Nishino2008.pdf')
(R/'qc/nishino_source_verification.json').write_text(json.dumps({'doi':'10.1254/jphs.08107FP','url':'https://www.jstage.jst.go.jp/article/jphs/107/3/107_08107FP/_pdf','local_pdf_reviewed':str(p),'pdf_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'table1_page353_visual_review':True,'all_80_drug_dose_time_cells_retained':True,'control_cells_separate':4,'reported_ED50_not_refitted':True,'flags':['triazolam 0.5mg/kg 90min 9/10 retained','reported ED50 time aggregation unspecified','parenthetical interval confidence level unspecified','reported diazepam ED50 3.11 is not corrected to match 60min 5/10 at 5mg/kg']},indent=2)+'\n')
print('Transcribed 80 drug-dose-time cells, 4 control cells, 5 reported ED50 values.')
