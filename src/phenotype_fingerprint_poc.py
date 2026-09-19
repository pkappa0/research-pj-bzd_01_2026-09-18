"""Fixed-X, within-study phenotype evidence PoC. No supervised learning or X selection."""
from pathlib import Path
import json, hashlib, shutil, itertools
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT.parent/'phenotype_poc'
DRUGS=['diazepam','alprazolam','triazolam','zolpidem','lorazepam','clonazepam','midazolam','temazepam','zopiclone','zaleplon']
AXES={'P1':'Motor coordination','P2':'Muscle relaxation','P3':'Sedation / locomotion','P4':'Hypnosis / LORR'}
KEYS=['study_id','assay_name','species','strain','sex','administration_route','dose','dose_unit','observation_time','outcome_name','vehicle','co_treatment','acute_chronic']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def csv(p,rows):pd.DataFrame(rows).to_csv(p,index=False,na_rep='NA')
def normalize(r):
    v=r['raw_value'];n=r['sample_size'];c=r['control_value']
    if pd.isna(v):return np.nan,'unavailable'
    if r['raw_unit']=='positive_mice_count' and pd.notna(n) and n>0:
        if not 0<=v<=n:raise ValueError('invalid failure count')
        return v/n,'positive_count / sample_size'
    if r['raw_unit']=='seconds_mean':return v,'raw duration seconds; no max scaling'
    if r['outcome_name'] in ['latency','locomotor_activity','grip_strength'] and pd.notna(c) and c>0:
        return (c-v)/c,'(control - drug) / control; no clipping'
    return np.nan,'no justified transform'
def rho_result(x,y):
    ok=np.isfinite(x)&np.isfinite(y);n=int(ok.sum())
    if n<3:return np.nan,n,'n_effective_lt_3'
    if len(set(x[ok]))<2 or len(set(y[ok]))<2:return np.nan,n,'constant_input'
    return float(spearmanr(x[ok],y[ok]).statistic),n,'descriptive_only'
def make_row(drug,axis,study,assay,**kw):
    row=dict(drug=drug,phenotype_axis=axis,assay_name=assay,species='mouse',strain='NA',sex='NA',administration_route='NA',dose=np.nan,dose_unit='mg/kg',observation_time='NA',outcome_name='qualitative_report',raw_value=np.nan,raw_unit='qualitative',control_value=np.nan,sample_size=np.nan,error_value=np.nan,error_type='NA',study_id=study,doi='NA',pmid='NA',source_url='NA',source_table_figure='abstract',extraction_method='manual primary-source paraphrase',evidence_quality='qualitative_primary',notes='',vehicle='NA',co_treatment='NA',acute_chronic='NA',original_record='NA',reported_direction='NA',result_available=True)
    row.update(kw);return row
def evidence(run):
    old=pd.read_csv(run/'raw/archived_mouse_phenotype.csv',keep_default_na=False);rows=[];audit=[]
    for i,r in old.iterrows():
        axis='P1' if r.assay=='rotarod' else 'P4' if 'righting' in r.assay else 'P3'
        num=lambda x:float(x) if str(x) else np.nan
        row=make_row(r.drug_id,axis,r.study_id,r.assay,strain=r.strain,administration_route=r.route,dose=num(r.dose_mg_per_kg),observation_time=str(r.time_min) if r.time_min else 'NA',outcome_name=r.endpoint,raw_value=num(r.raw_value),raw_unit=r.raw_unit,sample_size=num(r.n_tested),error_value=num(r.sd),error_type='archived_sd_unverified',doi=r.doi or 'NA',pmid=r.pmid,source_url=r.source_url,notes=r.condition_note,original_record=f'archived_mouse_phenotype.csv:{i+2}',reported_direction=r.effect_direction,extraction_method='archived extraction; source audit retained',evidence_quality='archived_quantitative' if r.raw_value else 'qualitative_primary')
        if r.study_id=='Nishino2008':
            row.update(strain='ICR',sex='male',vehicle='0.5% CMC sodium',co_treatment='none',acute_chronic='acute',control_value=0,source_table_figure='Table 1 p353; methods p349-350',evidence_quality='exact_table_checked')
        if r.study_id=='Tanaka2008':
            row.update(strain='ddY',sex='male',vehicle='saline' if r.drug_id=='zolpidem' else '0.5% PPG',acute_chronic='acute',source_table_figure='Table 1 p279' if axis=='P4' else 'Figure 1 / results p279')
            if axis=='P4':row.update(co_treatment='thiopental 20 mg/kg IV 10 min after drug',observation_time='duration after thiopental',administration_route='intraperitoneal',sample_size=7,control_value=113.1 if r.drug_id=='zolpidem' else 117.3,error_type='SEM',evidence_quality='exact_table_checked')
            else:row.update(co_treatment='none',observation_time='15/30/60 min; qualitative aggregate',sample_size=6)
        rows.append(row)
    audit += [dict(study='Nishino2008',field='strain; sex; control; vehicle',old='not_reported; missing',new='ICR; male; 0/10; 0.5% CMC sodium',source='Methods p349-350 and Table1 p353'),dict(study='Tanaka2008',field='error_type',old='sd',new='SEM',source='Table1 footnote p279'),dict(study='Tanaka2008',field='triazolam dose ambiguity',old='0.3 mg/kg archived',new='0.3 retained from Table1/results; methods paragraph says 3',source='Table1/results versus methods p278; unresolved source inconsistency')]
    for dose,v,err,n in [(3,146.3,41.6,7),(6,227.2,67.4,9)]:
        rows.append(make_row('zopiclone','P4','Tanaka2008','righting_reflex_with_thiopental',strain='ddY',sex='male',administration_route='intraperitoneal',dose=dose,observation_time='duration after thiopental',outcome_name='righting_reflex_duration',raw_value=v,raw_unit='seconds_mean',control_value=117.3,sample_size=n,error_value=err,error_type='SEM',doi='10.1254/jphs.FP0071991',pmid='18603831',source_url='https://www.jstage.jst.go.jp/article/jphs/107/3/107_FP0071991/_pdf',source_table_figure='Table1 p279',evidence_quality='exact_table_checked',vehicle='0.5% PPG',co_treatment='thiopental 20 mg/kg IV 10 min after drug',acute_chronic='acute',reported_direction='longer duration = more hypnosis',notes='Different doses and vehicles prohibit primary cross-drug matching.'))
    for drug,dose in [('zolpidem',3),('zopiclone',6),('triazolam',.3)]:
        for axis,assay in [('P1','rotarod'),('P2','traction')]:
            if axis=='P1' and drug!='zopiclone':continue
            rows.append(make_row(drug,axis,'Tanaka2008',assay,strain='ddY',sex='male',administration_route='intraperitoneal',dose=dose,observation_time='15/30/60 min; qualitative aggregate',sample_size=6,doi='10.1254/jphs.FP0071991',pmid='18603831',source_url='https://www.jstage.jst.go.jp/article/jphs/107/3/107_FP0071991/_pdf',source_table_figure='Figure1; results p279',vehicle='saline' if drug=='zolpidem' else '0.5% PPG',co_treatment='none',acute_chronic='acute',reported_direction='impairment reported' if drug=='triazolam' else 'no effect reported under tested conditions',notes='Qualitative; no numerical zero inferred.'))
    for d in ['zolpidem','zaleplon']:
        for axis,assay in [('P1','rotarod'),('P2','loaded_grid'),('P3','locomotor_activity')]:
            rows.append(make_row(d,axis,'Sanger1996',assay,pmid='8905326',doi='10.1016/0014-2999(96)00510-9',source_url='https://pubmed.ncbi.nlm.nih.gov/8905326/',reported_direction='deficit / suppression reported',notes='Abstract only. Loaded grid assigned P2 with mixed motor/strength construct warning; no exact dose/time extracted.'))
    for d in ['diazepam','lorazepam']:
        for assay in ['rotarod','beam_walking']:
            rows.append(make_row(d,'P1','Stanley2005',assay,pmid='15888506',doi='10.1177/0269881105051524',source_url='https://pubmed.ncbi.nlm.nih.gov/15888506/',reported_direction='impairment reported',notes='Receptor occupancy percentages are not phenotype scores and not normalized.'))
    rows.append(make_row('clonazepam','P3','Ono1976','open_field',administration_route='oral',pmid='986989',doi='10.1254/fpj.72.297',source_url='https://pubmed.ncbi.nlm.nih.gov/986989/',reported_direction='locomotion reduced',notes='Mouse-specific abstract result; rat result not transferred.'))
    rows.append(make_row('clonazepam','P1','Barnhill1990','rotarod',strain='CD-1',sex='male',administration_route='intraperitoneal',dose=2,observation_time='multiple times up to 14 h',pmid='2162948',source_url='https://pubmed.ncbi.nlm.nih.gov/2162948/',reported_direction='ataxia greater and longer in aged animals',notes='Age groups 6 weeks, 6 months, 1 and 2 years; no pooled numerical score.'))
    rows.append(make_row('temazepam','P3','Marshall1997','locomotor_activity',administration_route='oral',dose=10,pmid='9264062',doi='10.1016/s0091-3057(96)00132-3',source_url='https://pubmed.ncbi.nlm.nih.gov/9264062/',acute_chronic='daily for 7 days; acute challenge',reported_direction='tolerance to sedative effect reported',notes='Chronic regimen dose, not verified acute challenge dose; excluded from acute quantitative comparisons.'))
    for dose,count in [(90,6),(71.25,5),(58,5),(46.56,3),(37.38,2),(30,1)]:
        rows.append(make_row('midazolam','P4','Shi2024_PMC11559241','loss_of_righting_reflex',strain='C57BL/6',sex='male',administration_route='intraperitoneal',dose=dose,observation_time='within 10 min',outcome_name='LORR_count',raw_value=count,raw_unit='positive_mice_count',sample_size=6,doi='10.1186/s40001-024-02142-6',source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11559241/',source_table_figure='Table1; Animals; drug dose methods',evidence_quality='exact_html_table',acute_chronic='acute dose determination',co_treatment='none',reported_direction='higher fraction = more hypnosis',notes='Juvenile 20-22 day mice. Remimazolam comparator outside fixed 10 drugs; no within-panel Tier1 comparison.'))
    for i,r in enumerate(rows):r['evidence_id']=f'Y{i+1:04d}'
    csv(run/'tables/source_audit_corrections.csv',audit)
    csv(run/'tables/phenotype_evidence_raw.csv',rows)
    return pd.DataFrame(rows)

def analyze(run):
    cfg=json.loads((run/'config/phenotype_analysis_config.json').read_text());xrun=ROOT/cfg['X_run'];assert sha(xrun/cfg['X_matrix'])==cfg['X_sha256']
    e=evidence(run); e['value'],e['normalization_formula']=zip(*(normalize(r) for _,r in e.iterrows()))
    e['block_key']=e[KEYS].apply(lambda row: '|'.join('NA' if pd.isna(v) else str(v) for v in row),axis=1)
    sizes=e[e.value.notna()].groupby('block_key').drug.nunique()
    e['tier']=[1 if pd.notna(r.value) and sizes.get(r.block_key,0)>=2 and r.study_id=='Nishino2008' else 2 if pd.notna(r.value) else 4 for _,r in e.iterrows()]
    e['tier_reason']=np.where(e.tier==1,'same fully recorded study/assay/dose/time/route/vehicle; two drugs',np.where(e.tier==2,'numeric but unmatched dose/vehicle/age/study or single panel drug','qualitative only; numerical value NA'))
    e['primary_eligible']=(e.tier==1)&e.value.notna()
    csv(run/'tables/phenotype_normalized_within_study.csv',e.to_dict('records'))
    csv(run/'tables/phenotype_comparability_tiers.csv',e[['evidence_id','drug','phenotype_axis','block_key','tier','tier_reason','primary_eligible']].to_dict('records'))
    csv(run/'tables/phenotype_axis_mapping.csv',[dict(assay=a,phenotype_axis=p,direction='higher impairment score = stronger impairment',caveat='Assay is a proxy, not a pure physiological axis; loaded_grid mixes strength/coordination') for a,p in sorted(set(zip(e.assay_name,e.phenotype_axis)))])
    profile=[];missing=[]
    for d in DRUGS:
        row={'drug':d}
        for p in AXES:
            sub=e[(e.drug==d)&(e.phenotype_axis==p)]
            records=[dict(evidence_id=r.evidence_id,value=None if pd.isna(r.value) else r.value,raw_unit=r.raw_unit,source=r.source_url,study=r.study_id,tier=int(r.tier),block=r.block_key) for _,r in sub.iterrows()]
            row[p+'_records']=json.dumps(records) if records else 'NA';row[p+'_best_tier']=int(sub.tier.min()) if len(sub) else np.nan
            missing.append(dict(drug=d,phenotype_axis=p,evidence_rows=len(sub),numeric_rows=int(sub.value.notna().sum()),tier1_rows=int((sub.tier==1).sum()),status='NA' if not len(sub) else 'qualitative_only' if not sub.value.notna().any() else 'numeric_available'))
        profile.append(row)
    csv(run/'tables/phenotype_profile_10drug.csv',profile);csv(run/'tables/phenotype_missingness_summary.csv',missing)
    # Freeze Y before joining frozen X. Preserve all sources and all features.
    save(run/'config/PHENOTYPE_RESULTS_FROZEN.json',dict(frozen_at=datetime.now().astimezone().isoformat(),artifacts={str(p.relative_to(run)):sha(p) for p in sorted((run/'tables').glob('*.csv'))}))
    x=pd.read_csv(xrun/cfg['X_matrix']).set_index('drug_id').loc[DRUGS];features=list(x.columns);assert x.shape==(10,48)
    joined=e.merge(x.rename_axis('drug').reset_index(),on='drug',how='left',validate='many_to_one');csv(run/'tables/fingerprint_phenotype_joined.csv',joined.to_dict('records'))
    corr=[];ranks=[];pairs=[]
    primary=e[e.primary_eligible]
    for key,g in primary.groupby('block_key',sort=True):
        assert g.drug.nunique()==len(g)
        for f in features:
            a=x.loc[g.drug,f].to_numpy(float);b=g.value.to_numpy(float);rho,n,status=rho_result(a,b)
            corr.append(dict(block_key=key,phenotype_axis=g.phenotype_axis.iloc[0],feature=f,rho=rho,n_effective=n,status=status,evidence_ids=';'.join(g.evidence_id)))
            for (idx,r),rx,ry in zip(g.iterrows(),pd.Series(a).rank(),pd.Series(b).rank()):ranks.append(dict(block_key=key,drug=r.drug,feature=f,x_value=x.loc[r.drug,f],y_value=r.value,x_rank=rx,y_rank=ry,n_effective=n))
        for (_,r),(_,s) in itertools.combinations(g.iterrows(),2):
            for rep in ['alpha1','alpha2','E']:
                cols=features if rep=='E' else [f for f in features if f.startswith(rep+'|')];a=x.loc[r.drug,cols].to_numpy(float);b=x.loc[s.drug,cols].to_numpy(float)
                pairs.append(dict(block_key=key,phenotype_axis=r.phenotype_axis,drug1=r.drug,drug2=s.drug,representation=rep,structural_euclidean=float(np.linalg.norm(a-b)),structural_cosine=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b))),phenotype_abs_difference=abs(r.value-s.value),phenotype_similarity=1-abs(r.value-s.value),evidence_ids=r.evidence_id+';'+s.evidence_id,n_effective_drugs=2,pair_independence='same drug pair repeated across conditions; not independent'))
    for axis in AXES:
        if not any(r['phenotype_axis']==axis for r in corr):corr += [dict(block_key='no_Tier1_numeric_context',phenotype_axis=axis,feature=f,rho=np.nan,n_effective=0,status='no_comparable_data',evidence_ids='NA') for f in features]
    csv(run/'tables/feature_phenotype_spearman.csv',corr);csv(run/'tables/feature_phenotype_ranks.csv',ranks);csv(run/'tables/structural_vs_phenotype_similarity.csv',pairs)
    dose=[]
    # Only within each drug / study / assay / time / vehicle, over observed doses.
    keys=[k for k in KEYS if k not in ['dose','dose_unit']]+['drug','phenotype_axis']
    for key,g in e[e.value.notna()].groupby(keys,dropna=False):
        g=g.sort_values('dose');vals=g.value.to_numpy(float);ds=g.dose.to_numpy(float)
        dose.append(dict(zip(keys,key),tested_doses=json.dumps(ds.tolist()),dose_unit='mg/kg',maximum_observed_effect=float(max(vals)),auc_observed_range=float(np.trapezoid(vals,ds)) if len(g)>1 else np.nan,auc_unit=('fraction' if g.raw_unit.iloc[0]=='positive_mice_count' else 'seconds')+' * mg/kg',lowest_effective_dose=np.nan,ED50=np.nan,notes='No significance-based threshold inferred. AUC descriptive within this drug/range only; different ranges not compared. No dose-only potency ranking.'))
    dose += [dict(drug=d,study_id='Nishino2008',phenotype_axis='P1',ED50=v,ED50_interval=ci,dose_unit='mg/kg',notes='Source-reported probit estimate, not refit; rotarod context only, no exposure adjustment') for d,v,ci in [('diazepam',3.11,'2.53-3.81'),('triazolam',1.25,'1.00-1.49')]]
    dose += [dict(drug='midazolam',study_id='Shi2024_PMC11559241',phenotype_axis='P4',ED50=46,ED50_interval='32.809-56.841',dose_unit='mg/kg',notes='Source-reported ED50; juvenile LORR, no cross-study potency comparison')]
    csv(run/'tables/within_study_dose_response.csv',dose)
    figures(run,x,e,pd.DataFrame(corr),pd.DataFrame(pairs))
    baseline=json.loads((WORK/'baseline.json').read_text());assert all(sha(ROOT/p)==h for p,h in baseline.items())
    assert sha(xrun/cfg['X_matrix'])==cfg['X_sha256']
    frozen=json.loads((xrun/'STRUCTURAL_RESULTS_FROZEN.json').read_text());assert all(sha(xrun/p)==h for p,h in frozen['artifacts'].items())
    qc=dict(status='PASS_WITH_EVIDENCE_LIMITATIONS',evidence_rows=len(e),drugs_with_any_evidence=e.drug.nunique(),drugs_with_numeric_evidence=e[e.value.notna()].drug.nunique(),axis_coverage={p:dict(any=e[e.phenotype_axis==p].drug.nunique(),numeric=e[(e.phenotype_axis==p)&e.value.notna()].drug.nunique(),tier1=e[(e.phenotype_axis==p)&(e.tier==1)].drug.nunique()) for p in AXES},tier1_drugs=primary.drug.nunique(),tier1_contexts=primary.block_key.nunique(),unique_Tier1_drug_pairs=len(set(tuple(sorted((r['drug1'],r['drug2']))) for r in pairs)),estimable_spearman=sum(pd.notna(r['rho']) for r in corr),all_48_X_features_retained=True,X_hash_unchanged=True,old_files_hash_unchanged=len(baseline),no_docking_or_clustering=True,no_supervised_ML=True,no_missing_zero_imputation=True,ML_rationale='indeterminate')
    save(run/'reports/QC_REPORT.json',qc);report(run,qc)
    shutil.copy2(__file__,run/'code/phenotype_fingerprint_poc.py')
    save(run/'run_manifest.json',dict(completed_at=datetime.now().astimezone().isoformat(),artifacts={str(p.relative_to(run)):sha(p) for p in sorted(run.rglob('*')) if p.is_file() and p.name!='run_manifest.json'}))
    print(json.dumps(qc,indent=2))

def figures(run,x,e,corr,pairs):
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    def out(fig,name):
        fig.savefig(run/'figures'/f'{name}.png',dpi=145,bbox_inches='tight');fig.savefig(run/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,5));ax.axis('off');ax.text(.05,.5,'Unsteadiness\n(multifactorial phenotype)',ha='left',va='center',fontsize=16,bbox=dict(boxstyle='round',fc='#e8edf4'))
    for i,(p,label) in enumerate(AXES.items()):
        yy=.86-i*.23;ax.annotate(p+'  '+label,xy=(.36,.5),xytext=(.55,yy),ha='left',va='center',arrowprops=dict(arrowstyle='<-',color='#596d8c'),bbox=dict(boxstyle='round',fc='#edf4ef'),fontsize=13)
    fig.suptitle('Separate measured axes; no composite score or causal decomposition');out(fig,'phenotype_decomposition')
    tier=np.full((10,4),np.nan);labels={}
    for i,d in enumerate(DRUGS):
        for j,p in enumerate(AXES):
            sub=e[(e.drug==d)&(e.phenotype_axis==p)];tier[i,j]=sub.tier.min() if len(sub) else np.nan;labels[i,j]='NA' if not len(sub) else f'T{int(tier[i,j])}\n{sub.value.notna().sum()} numeric / {len(sub)} records'
    fig,ax=plt.subplots(figsize=(12,7));cmap=plt.get_cmap('viridis_r').copy();cmap.set_bad('#ddd');ax.imshow(tier,vmin=1,vmax=4,cmap=cmap,aspect='auto');ax.set_xticks(range(4),[p+' '+a for p,a in AXES.items()],fontsize=9);ax.set_yticks(range(10),DRUGS)
    for (i,j),txt in labels.items():ax.text(j,i,txt,ha='center',va='center',fontsize=9,color='white' if tier[i,j]==4 else 'black')
    ax.set_title('Evidence coverage, not effect size | T4 qualitative; gray NA | T1 fully matched');out(fig,'phenotype_evidence_matrix')
    fig,ax=plt.subplots(figsize=(10,17));ax.imshow(np.full((48,4),np.nan),vmin=-1,vmax=1,cmap=cmap,aspect='auto');ax.set_yticks(range(48),[f.replace('|any_contact_frequency','').replace('BZD_GAMMA2_','gamma2_').replace('BZD_SITE_','alpha_') for f in x.columns],fontsize=8);ax.set_xticks(range(4),list(AXES))
    for i in range(48):
        for j,p in enumerate(AXES):ax.text(j,i,'NA\nn='+str(int(corr[corr.phenotype_axis==p].n_effective.max())),ha='center',va='center',fontsize=7)
    ax.set_title('All 48 frozen features retained\nSpearman not estimable: no Tier1 context has >=3 drugs');out(fig,'fingerprint_phenotype_correlation_heatmap')
    g=e[(e.study_id=='Nishino2008')&(e.dose==2)&(e.observation_time=='60.0')].sort_values('value');fig,(ax,ay)=plt.subplots(1,2,figsize=(18,4),gridspec_kw={'width_ratios':[8,1]});ax.imshow(x.loc[g.drug],vmin=0,vmax=1,cmap='viridis',aspect='auto');ax.set_yticks(range(len(g)),g.drug);ax.set_xticks(range(48),[f.split('|')[0]+':'+f.split('|')[2]+f.split('|')[1].split('_')[-1] for f in x.columns],rotation=90,fontsize=6);ax.axvline(20.5,color='white',lw=2);ax.text(.5,1.02,'Contact frequency: purple = 0, yellow = 1',transform=ax.transAxes,ha='center',fontsize=9);ay.barh(range(len(g)),g.value,color='#8c5849');ay.invert_yaxis();ay.set_xlim(0,1);ay.set_yticks([]);ay.set_xlabel('Failure fraction');fig.suptitle('Nishino2008 | oral 2 mg/kg, 60 min | n_drugs=2 | all frozen alpha1 (21) + alpha2 (27) features');out(fig,'phenotype_sorted_fingerprint')
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    for ax,rep in zip(axs,['alpha1','alpha2','E']):
        p=pairs[pairs.representation==rep];ax.scatter(p.structural_cosine,p.phenotype_similarity,c=np.arange(len(p)),cmap='viridis',s=70);ax.set_xlim(0,1.02);ax.set_ylim(0,1.02);ax.set_xlabel('Fixed fingerprint cosine similarity');ax.set_ylabel('1 - absolute failure-fraction difference');ax.set_title(rep+' | one unique drug pair')
    fig.suptitle('Eight dose/time contexts repeat diazepam-triazolam; no association estimate or fit');out(fig,'structural_vs_phenotype_similarity')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for ax,feature in zip(axs,[x.columns[1],x.columns[22]]):
        for _,r in g.iterrows():ax.scatter(x.loc[r.drug,feature],r.value);ax.annotate(r.drug,(x.loc[r.drug,feature],r.value),xytext=(4,4),textcoords='offset points',fontsize=8)
        ax.set_xlabel(feature.replace('|any_contact_frequency',''));ax.set_ylabel('P1 failure fraction');ax.set_xlim(-.05,1.05);ax.set_ylim(-.05,1.05);ax.set_title('Prespecified TYR58 example; n=2, rho=NA')
    fig.tight_layout();out(fig,'feature_phenotype_scatter')

def report(run,q):
    text='''# Fixed fingerprint × phenotype decomposition PoC

既存10剤の構造fingerprintを変更せず、Y側をP1 motor coordination、P2 muscle relaxation、P3 sedation/locomotion、P4 hypnosisに分離した。各assayは複数の生理機能に影響されるproxyであり、4つの純粋な独立機序を測定したとは解釈しない。

## A. Phenotype feasibility
''' + f"対象10剤中、定性を含む文献証拠は{q['drugs_with_any_evidence']}剤、数値は{q['drugs_with_numeric_evidence']}剤。raw evidence {q['evidence_rows']}行。coverageは次表。\n\n| Axis | Any evidence | Numeric | Tier1 drugs |\n|---|---:|---:|---:|\n" + '\n'.join(f"|{p}|{r['any']}|{r['numeric']}|{r['tier1']}|" for p,r in q['axis_coverage'].items()) + '''

Tier1はNishino2008のdiazepam/triazolam、同一用量2/5 mg/kg × 15/30/60/90分の8 contextのみ。薬剤数は2、unique drug pairは1。時間・用量をreplicate drugsにしない。Nishino表のstrain/sex/controlを原論文から補足。Tanaka2008は誤差表記をSEMとして新run内で訂正し、旧ファイルは不変。Table1のtriazolam0.3 mg/kgとmethods本文3 mg/kgの不一致も明示した。

TanakaのP4はthiopental併用のdurationで、単独薬hypnosisとは別。用量・saline/PPGが異なりTier1にしない。Midazolamは幼若雄C57BL/6のLORR countであり、他studyへ転用しない。Temazepamは反復投与後の耐性の定性記録で、急性効果量ではない。T4「効果なし」は0ではなく数値NA。phenotype_profileは各cellの全value/source/tier/contextをJSONで保持し、study横断の代表値を作らない。

正規化は失敗数/n。durationは秒のまま保持。controlとの差によるreductionの実装はあるが、対応するcontrolと実数値がない場合は計算しない。元表の全dose-responseを保持し、同一drug・study・time・assay内のobserved-range AUCとmaximumを記述。試験用量域が違うAUCを比較しない。最低有効用量はsignificanceを再構成せずNA。原論文ED50は出典付きで別欄に保持し、mg/kgから薬効順位を作らない。

## B. Structural correspondence
**indeterminate**。全48個の固定X featureに対してcontextごとのSpearman表を保存したが、Tier1はn_effective=2であり、rhoはNAとした。2点なら非tie時に必ず±1となるため、高相関候補として報告しない。P2/P3/P4にはTier1数値contextがなくn=0。featureの選別や削除はしていない。

P1ではtriazolamのfailure fractionがdiazepamより8条件中7条件で高く1条件で同値。これは同じ1対についての方向の記述に限る。rank表・全feature heatmap・事前指定TYR58のscatterを保存。特定残基やfeature combinationがphenotypeを説明するとの候補確定には足りない。「最も強く対応するaxis」は判断不能で、P1は最も比較可能な数値があるaxisにすぎない。

## C. Multi-receptor relevance
**indeterminate**。同じ1対の構造類似度は各block内で固定値なのにphenotype差はdose/timeで変わる。Figure5の8点は独立8drug-pairsではない。α1、α2、連結Eの3パネルを並べたが、どれがphenotypeをよりよく説明するかの比較統計は成立しない。連結の有用性を示したとは結論しない。

## D. ML rationale
**indeterminate**。4axisの情報整理は実行可能だが、比較可能なYを備えたdrug数が足りない。教師ありモデルは実行していない。次の反証可能な仮説は、同一assay/time/routeと測定された脳内曝露を揃えた独立薬剤群で、固定Eの薬剤間距離がP1/P2/P3/P4それぞれの差と関連し、single-blockや単純化学記述子を超える情報を持つか、である。再現しなければ追加情報仮説は支持されない。

## Limitations / source audit
最大の制約は定性coverageと比較可能な数値coverageの乖離：10剤に文献証拠があってもTier1は2剤のみ。PK、年齢、溶媒、併用薬、投与歴、species、endpointを跨いだ混合はしていない。検索は記録した公開情報の範囲でありsystematic reviewではない。全文未取得の論文では数値抽出を推測しない。Mouse以外・4axis外の候補は別表で保持。既存のP1/P2というevidence labelは今回のphenotype axisやcomparability tierとは無関係に再評価した。

構造側にはβ3構造、α2 provisional local construct、10剤のみ、単一conformer・固定環・未指定立体・pH microstate未列挙が残る。動物phenotypeからヒトふらつき・副作用を説明したとはいえない。

## Sources
- [Nishino2008](https://doi.org/10.1254/jphs.08107FP): Table1、methods、reported ED50。
- [Tanaka2008](https://doi.org/10.1254/jphs.FP0071991): Table1、Figure1とmethods/results。
- [Sanger1996](https://pubmed.ncbi.nlm.nih.gov/8905326/)、[Stanley2005](https://pubmed.ncbi.nlm.nih.gov/15888506/)、[Ono1976](https://pubmed.ncbi.nlm.nih.gov/986989/)、[Barnhill1990](https://pubmed.ncbi.nlm.nih.gov/2162948/)、[Marshall1997](https://pubmed.ncbi.nlm.nih.gov/9264062/): abstract-level qualitative evidence。
- [Shi2024](https://doi.org/10.1186/s40001-024-02142-6): Table1、juvenile methods。既存Bayley/Bourin/Henauer/Lopezも元行とsource IDを保持。

![Evidence matrix](../figures/phenotype_evidence_matrix.png)
![Phenotype-sorted fixed fingerprint](../figures/phenotype_sorted_fingerprint.png)
![Correspondence availability](../figures/fingerprint_phenotype_correlation_heatmap.png)
'''
    (run/'reports/FINAL_REPORT.md').write_text(text)

if __name__=='__main__':analyze(Path((WORK/'active_run.txt').read_text()))
