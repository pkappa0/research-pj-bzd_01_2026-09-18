"""Immutable, descriptive X-y prototype. No fitting or cross-validation.

Read ANALYSIS_CONTRACT.md before interpretation. Primary comparisons fix receptor.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime
from fractions import Fraction
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from . import receptor_fixed_drug_comparison as rf

ROOT = rf.ROOT
DRUGS, BLOCKS = rf.DRUGS, rf.BLOCKS


def csv_write(path, rows):
    """Explicit NA on disk, rather than an ambiguous empty numeric field."""
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader()
        for r in rows:
            w.writerow({k: 'NA' if r.get(k) is None else json.dumps(r[k], ensure_ascii=False, sort_keys=True) if isinstance(r[k], (list, dict)) else r[k] for k in keys})


def display_keep(d, a, config):
    if d is None or a is None:
        return True
    return not (d < Fraction(config['frequency_threshold']) and a < Fraction(config['frequency_threshold'])
                and abs(a-d) <= Fraction(config['maximum_absolute_difference']))


def frequency_fraction(row, metric):
    count = row['n_interacting_poses'] if metric == 'frequency' else row[metric.replace('_frequency','_n_interacting_poses')]
    value = row[metric]
    if value is None:
        return None
    if count is None or not row['n_poses'] or not 0 <= count <= row['n_poses']:
        raise ValueError('Invalid audited frequency counts')
    result = Fraction(count, row['n_poses'])
    if abs(float(result)-value) > 1e-12:
        raise ValueError('Frequency does not match contacting-pose counts')
    return result


def prepare(config):
    for key in ['ifp_run', 'pharmacology_run', 'reference_run']:
        rf.verify_manifest(ROOT/config[key])
    inputs = rf.read_jsonl(ROOT/config['ifp_run']/'processed/diazepam_alprazolam_ifp.jsonl')
    activities = rf.read_jsonl(ROOT/config['pharmacology_run']/'processed/bzd_activity_classification.jsonl')
    structural = rf.structural_comparison(inputs, config['compound_ids'])
    pharm, selected, unpaired = rf.pharmacology_comparison(activities, config['compound_ids'])
    pharm.sort(key=lambda r: r['receptor_block'])
    raw, schema = rf.multi_receptor_matrix(structural)
    index = {(r['receptor_subtype'], r['drug_id'], r['feature_id']):r for r in inputs}
    audit = []
    for row, col in zip(structural, schema):
        base = row['feature'].rsplit('__',1)[0]
        ds, als = [index[(row['receptor_block'], d, base)] for d in DRUGS]
        metric = row['metric'] if row['feature_type'] != 'geometry' else 'frequency'
        d, a = frequency_fraction(ds,metric), frequency_fraction(als,metric)
        keep = display_keep(d,a,config)
        audit.append(dict(col, retained_for_presentation=keep,
            selection_basis='parent_pi_frequency' if row['feature_type']=='geometry' else 'own_frequency',
            diazepam_exact_frequency=str(d), alprazolam_exact_frequency=str(a),
            exact_absolute_difference=str(abs(a-d)) if d is not None and a is not None else None,
            reason='retain_missing_evidence' if d is None or a is None else 'retained_by_config' if keep else 'both_low_and_small_difference'))
    keys = [r['feature_key'] for r in audit if r['retained_for_presentation']]
    reduced = [{k:v for k,v in r.items() if k=='drug_id' or k in keys} for r in raw]
    targets, yschema = [], []
    for p in pharm:
        # Every observed endpoint gets its own context-qualified column; none pooled.
        key = 'y|'+p['receptor_composition']+'|'+p['species']+'|'+p['endpoint_class']+'|'+p['endpoint']+'|'+p['unit']+'|'+p['diazepam_assay_chembl_id']
        yschema.append(dict(column=key, receptor_block=p['receptor_block'], receptor_composition=p['receptor_composition'], species=p['species'], endpoint_class=p['endpoint_class'], endpoint=p['endpoint'], unit=p['unit'], relation=p['relation'], assay_id=p['diazepam_assay_chembl_id'], document_id=p['document'], availability='observed_same_recorded_assay_pair'))
        for drug in DRUGS:
            s=next(s for s in selected if s['activity_id']==p[drug+'_activity_id'])
            targets.append(dict(drug_id=drug, molecule_chembl_id=s['molecule_chembl_id'], y_column=key,
                receptor_block=p['receptor_block'], receptor_composition=s['receptor_composition'], species=s['species'],
                endpoint_class=s['endpoint_class'], endpoint=s['standard_type'], value=s['standard_value'], unit=s['standard_units'], relation=s['standard_relation'],
                assay_id=s['assay_chembl_id'], document_id=s['document_chembl_id'], source_activity_id=s['activity_id'], target_id=s['target_chembl_id'],
                assay_description=s['assay_detail_description'], assay_parameters=s['assay_detail_assay_parameters'],
                source_run=config['pharmacology_run'], acquisition_run=s['classification_source_run'], source_file=s['source_file'], source_json_pointer=s['source_json_pointer'], source_sha256=s['source_sha256'],
                pharmacology_context_match_status=p['pharmacology_context_match_status'], structure_pharmacology_context_match_status=p['structure_pharmacology_context_match_status']))
    for block in BLOCKS:
        for family, endpoint in [('binding','Kd'),('binding','other_affinity'),('functional_potency','EC50'),('functional_potency','IC50'),('functional_efficacy','Emax'),('functional_efficacy','response')]:
            yschema.append(dict(column=f'y|{block}_beta2_gamma2|human|{family}|{endpoint}|unavailable', receptor_block=block, receptor_composition=block+'_beta2_gamma2', species='human', endpoint_class=family, endpoint=endpoint, unit=None, relation=None, assay_id=None, document_id=None, availability='no_comparable_two_drug_pair_in_selected_source; context_is_requested_not_observed'))
    dataset, mappings = [], []
    for rawrow in raw:
        drug=rawrow['drug_id']; record=dict(rawrow, molecule_chembl_id=config['compound_ids'][drug])
        for y in yschema:
            found=[t for t in targets if t['drug_id']==drug and t['y_column']==y['column']]
            t=found[0] if found else None
            record[y['column']]=t['value'] if t else None
            status=t['structure_pharmacology_context_match_status'] if t else 'not_assessed_missing_activity'
            record['qc|'+y['column']]=status
            mappings.append(dict(drug_id=drug, molecule_chembl_id=config['compound_ids'][drug], x_block=y['receptor_block'], x_columns=[s['feature_key'] for s in schema if s['receptor_block']==y['receptor_block']],
                y_column=y['column'], value=record[y['column']], unit=y['unit'], endpoint=y['endpoint'], endpoint_class=y['endpoint_class'],
                structural_context='alpha1_beta3_gamma2L_full_pentamer' if y['receptor_block']=='alpha1' else 'provisional_beta3_alpha2_gamma2_local_CDE',
                pharmacology_context=y['receptor_composition'], species=y['species'], context_match_status=status,
                context_reason='known_beta2_vs_beta3_mismatch; alpha2 additionally local construct' if t else 'No observed y; planned context only',
                pharmacology_context_match_status=t['pharmacology_context_match_status'] if t else None,
                source_activity_id=t['source_activity_id'] if t else None, assay_id=y['assay_id'], document_id=y['document_id']))
        dataset.append(record)
    # Supplemental model 0-2 candidates: pose unions, NOT a sum of type frequencies.
    baseline = [{ 'drug_id':d } for d in DRUGS]
    docking=list(csv.DictReader((ROOT/config['docking_source']).open()))
    for block in BLOCKS:
        sites=list(dict.fromkeys(r['site_label'] for r in structural if r['receptor_block']==block))
        for drug, record in zip(DRUGS,baseline):
            rows=[r for r in inputs if r['receptor_subtype']==block and r['drug_id']==drug]
            for site in sites:
                members=[r for r in rows if r['site_label']==site]
                if {r['n_poses'] for r in members}!={9}:raise ValueError('Unequal denominator')
                poses=set(p for r in members for p in r['pose_ids'])
                record[f'{block}|{site}__any_selected_contact__binary']=int(bool(poses))
                record[f'{block}|{site}__any_selected_contact__frequency']=len(poses)/9
            scores=[float(r['vina_score_kcal_mol']) for r in docking if r['drug_id']==drug and r['receptor_id']==rows[0]['receptor_id']]
            if len(scores)!=9:raise ValueError('Expected nine archived docking scores')
            record[f'{block}|vina_score_min_kcal_mol']=min(scores)
    manifest=[]
    for model in range(6):
        for block in (BLOCKS if model<5 else ('alpha1|alpha2',)):
            if model==0:columns=[block+'|vina_score_min_kcal_mol'];source='baseline_features.csv'
            elif model in [1,2]:columns=[k for k in baseline[0] if k.startswith(block+'|') and k.endswith('__binary' if model==1 else '__frequency')];source='baseline_features.csv'
            else:columns=[s['feature_key'] for s in schema if (model==5 or s['receptor_block']==block) and (model>=4 or s['feature_type']!='geometry')];source='raw_fingerprint.csv'
            manifest.append(dict(model=model, receptor_scope=block, feature_table=source, columns=columns, feature_count=len(columns),
                representation=['docking minimum archived score','binary any selected residue contact','type-collapsed pose-union frequency','type-resolved frequency plus P/T','type-resolved frequency plus P/T plus geometry','concatenated alpha1 and alpha2 model4'][model],
                target_policy='same endpoint, assay/context, drug split and eligible drug cohort across models; evaluate each y separately',
                future_evaluation='drug-grouped held-out evaluation after n expansion; preprocessing fit within training folds; compare model5 against best prespecified single-block baseline',
                state='definition_and_features_only; no_training_no_CV_no_importance',
                limitations='Model1/2 cover six named sites and selected interaction types only; not all receptor residues. Docking poses not independent drugs. Geometry missingness retained.'))
    return dict(inputs=inputs, structural=structural, pharm=pharm, selected=selected, unpaired=unpaired, raw=raw, schema=schema, audit=audit, reduced=reduced, targets=targets, yschema=yschema, dataset=dataset, mappings=mappings, baseline=baseline, manifest=manifest)


def label(s):
    return s.replace('alpha1|','α1 ').replace('alpha2|','α2 ').replace('gamma2_','γ2 ').replace('alpha_','α ').replace('__hydrophobic_interaction__frequency',' hyd').replace('__hydrogen_bond__frequency',' H-bond').replace('__halogen_bond__frequency',' halogen').replace('__pi_stack__frequency',' π').replace('__pi_stack__type_P_frequency',' π/P').replace('__pi_stack__type_T_frequency',' π/T')


def figures(path, data):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    def save(fig,name):
        for ext in ['png','pdf']:fig.savefig(path/(name+'.'+ext),dpi=170,bbox_inches='tight')
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(15,7));ax.axis('off');ax.set(xlim=(0,15),ylim=(0,7))
    ax.text(7.5,6.65,'Structural Fingerprint → Pharmacological Activity Profile',ha='center',fontsize=19,weight='bold')
    ax.text(1,5.9,'X = residue interaction fingerprint',fontsize=13,color='#235779')
    ax.text(9,5.9,'y = pharmacological activity database',fontsize=13,color='#805216')
    for i,drug in enumerate(DRUGS):
        y=4.7-2.1*i
        ax.text(.1,y,drug.title(),weight='bold',fontsize=13)
        for x,txt,color in [(2,'α1 block\n27 features','#dfedf7'),(4.7,'α2 block*\n27 features','#e0efdf')]:
            ax.text(x,y,txt,bbox=dict(boxstyle='round,pad=.6',fc=color,ec='#667'),va='center',fontsize=12)
        ax.annotate('',xy=(8.8,y),xytext=(7.2,y),arrowprops=dict(arrowstyle='->',lw=2))
        vals=[next(t['value'] for t in data['targets'] if t['drug_id']==drug and t['receptor_block']==b) for b in BLOCKS]
        ax.text(9,y,f'α1β2γ2 Ki: {float(vals[0]):g} nM\nα2β2γ2 Ki: {float(vals[1]):g} nM\nFunctional potency / efficacy: NA',va='center',fontsize=12,bbox=dict(boxstyle='round,pad=.6',fc='#fff0d9',ec='#997'))
    ax.text(7.5,.7,'Data linkage only — no model trained (n = 2 drugs)\nPrimary axis: different drugs within each receptor block. β3 structure ↔ β2 pharmacology: context_mismatch.\n*α2 is a provisional local construct. Geometry NA is retained; X blocks are concatenated, never averaged.',ha='center',fontsize=10)
    save(fig,'structural_fingerprint_to_activity_profile')
    # Matrix panels use distinct scales, including distinct geometry units.
    groups=[('X: frequency (0–1)',[s for s in data['schema'] if s['feature_type']!='geometry']),('X: centroid distance (Å)',[s for s in data['schema'] if s['feature_key'].endswith('centdist_mean')]),('X: angle (degrees)',[s for s in data['schema'] if s['feature_key'].endswith('angle_mean')]),('X: offset (Å)',[s for s in data['schema'] if s['feature_key'].endswith('offset_mean')])]
    fig,axes=plt.subplots(5,1,figsize=(18,16),layout='constrained')
    for ax,(title,cols) in zip(axes,groups):
        arr=np.array([[np.nan if r[s['feature_key']] is None else float(r[s['feature_key']]) for s in cols] for r in data['raw']])
        im=ax.imshow(np.ma.masked_invalid(arr),aspect='auto',cmap=plt.get_cmap('YlGnBu').with_extremes(bad='#dddddd'),vmin=0,vmax=1 if 'frequency' in title else None)
        ax.set(title=title,yticks=[0,1],yticklabels=DRUGS,xticks=range(len(cols)),xticklabels=[label(s['feature_key']).replace('__pi_stack__',' π ') for s in cols])
        ax.tick_params(axis='x',labelrotation=65,labelsize=8);fig.colorbar(im,ax=ax,shrink=.65)
        for i,j in zip(*np.where(np.isnan(arr))):ax.text(j,i,'NA',ha='center',va='center',fontsize=8)
    ax=axes[-1]; cols=data['yschema']; arr=np.array([[float(r[y['column']]) if r[y['column']] is not None else np.nan for y in cols] for r in data['dataset']])
    im=ax.imshow(np.ma.masked_invalid(arr),aspect='auto',cmap=plt.get_cmap('Oranges').with_extremes(bad='#ddd'),vmin=0)
    ax.set(title='y: observed Ki (nM); all other endpoint slots NA, not on the Ki scale',yticks=[0,1],yticklabels=DRUGS,xticks=range(len(cols)),xticklabels=[y['receptor_block']+' '+y['endpoint'] for y in cols]);ax.tick_params(axis='x',labelrotation=45)
    for i in range(2):
        for j in range(len(cols)):ax.text(j,i,'NA' if np.isnan(arr[i,j]) else f'{arr[i,j]:g}',ha='center',va='center')
    fig.colorbar(im,ax=ax,label='Observed Ki only (nM)',shrink=.6)
    fig.suptitle('Two-drug dataset prototype · block-qualified X, endpoint-qualified y\nGray = NA; geometry conditional on observed contacts; β2/β3 context mismatch',fontsize=16)
    save(fig,'ml_dataset_matrix')
    freq=[s for s in data['audit'] if s['retained_for_presentation'] and s['feature_type']!='geometry']
    fig,axes=plt.subplots(1,3,figsize=(17,5),gridspec_kw={'width_ratios':[3,1,1]},layout='constrained')
    arr=np.array([[r[s['feature_key']] for s in freq] for r in data['raw']])
    im=axes[0].imshow(arr,vmin=0,vmax=1,aspect='auto',cmap='YlGnBu');axes[0].set(yticks=[0,1],yticklabels=DRUGS,xticks=range(len(freq)),xticklabels=[label(s['feature_key']) for s in freq],title='Mechanically retained frequency features')
    axes[0].tick_params(axis='x',labelrotation=60)
    for i in range(2):
        for j in range(len(freq)):axes[0].text(j,i,f'{arr[i,j]:.2f}',ha='center',va='center',color='white' if arr[i,j]>.6 else 'black')
    fig.colorbar(im,ax=axes[0],label='Pose frequency',shrink=.65)
    for ax,p in zip(axes[1:],data['pharm']):
        v=[float(p[d+'_value']) for d in DRUGS];ax.bar([0,1],v,color=['#315f83','#a35c28']);ax.set(xticks=[0,1],xticklabels=DRUGS,ylabel='Ki (nM)',title=p['receptor_composition'].replace('_','\n'),ylim=(0,max(v)*1.3));ax.tick_params(axis='x',labelrotation=35)
        for i,val in enumerate(v):ax.text(i,val+max(v)*.03,f'{val:g}',ha='center')
    fig.suptitle('Presentation filter only · retained does not mean significant or drug-discriminating\nSmaller Ki = higher reported affinity; structural/pharmacology contexts mismatch',fontsize=13)
    save(fig,'reduced_fingerprint_activity_summary')


def run(config_path, tests_xml):
    tracked=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
    before={f:rf.sha((ROOT/f).read_bytes()) for f in tracked if f}
    config=json.loads(config_path.read_text()); data=prepare(config)
    test_summary=rf.tests_summary(tests_xml)
    path=ROOT/'runs'/(datetime.now().astimezone().strftime('%Y%m%d_%H%M')+'_structure_activity_ml_prototype')
    path.mkdir(exist_ok=False)
    for name in ['tables','figures','reports','config','provenance','code']: (path/name).mkdir()
    tables={'raw_fingerprint':data['raw'],'reduced_fingerprint':data['reduced'],'pharmacology_targets':data['targets'],
        'ml_dataset_prototype':data['dataset'],'receptor_fixed_drug_differences':data['structural'],
        'receptor_fixed_pharmacology_differences':data['pharm'],'structure_activity_mapping':data['mappings'],
        'feature_set_manifest':data['manifest'],'feature_schema':data['schema'],'target_schema':data['yschema'],
        'presentation_filter_audit':data['audit'],'baseline_features':data['baseline']}
    major=[s['feature_key'] for s in data['audit'] if s['retained_for_presentation'] and s['feature_type']!='geometry']
    tables['presentation_activity_summary']=[{k:v for k,v in r.items() if k in ['drug_id']+major+[y['column'] for y in data['yschema'] if y['availability'].startswith('observed')]} for r in data['dataset']]
    tables['missingness_mask']=[{k:v is None for k,v in r.items() if not k.startswith('qc|')} | {'drug_id':r['drug_id']} for r in data['dataset']]
    for name,rows in tables.items():csv_write(path/'tables'/(name+'.csv'),rows)
    for name,rows in [('source_ifp',data['inputs']),('selected_activity_full',data['selected']),('eligible_unpaired_activities',data['unpaired'])]:
        (path/'provenance'/(name+'.jsonl')).write_text(''.join(rf.encode(r)+'\n' for r in rows))
    sources=[ROOT/config['ifp_run']/'processed/diazepam_alprazolam_ifp.jsonl',ROOT/config['pharmacology_run']/'processed/bzd_activity_classification.jsonl',ROOT/config['docking_source'],ROOT/'ANALYSIS_CONTRACT.md',ROOT/config['reference_run']/'multi_receptor_drug_fingerprint.csv']
    rf.write_json(path/'provenance/source_inventory.json',[dict(path=str(p.relative_to(ROOT)),sha256=rf.sha(p.read_bytes())) for p in sources])
    shutil.copy2(ROOT/'ANALYSIS_CONTRACT.md',path/'provenance/ANALYSIS_CONTRACT.md')
    shutil.copy2(config_path,path/'config/fingerprint_reduction_rules.json');shutil.copy2(tests_xml,path/'reports/tests.xml')
    for name in ['structure_activity_ml_prototype.py','receptor_fixed_drug_comparison.py','bzd_pharmacology_filter.py','run_manager.py']:
        shutil.copy2(ROOT/'src'/name,path/'code'/name)
    shutil.copy2(ROOT/'tests/test_structure_activity_ml_prototype.py',path/'code/test_structure_activity_ml_prototype.py')
    figures(path/'figures',data)
    reference=list(csv.DictReader(sources[-1].open()))
    raw_matches=all((r[k] is None and ref[k]=='') or (r[k] is not None and float(r[k])==float(ref[k])) for r,ref in zip(data['raw'],reference) for k in r if k!='drug_id')
    checks=dict(preexisting_tracked_files_unchanged=all(rf.sha((ROOT/f).read_bytes())==digest for f,digest in before.items()), raw_ifp_matches_previous_run=raw_matches, activity_values_reextracted=len(data['targets'])==4 and all(t['value']==next(s['standard_value'] for s in data['selected'] if s['activity_id']==t['source_activity_id']) for t in data['targets']),
        endpoint_columns_separate=len({y['column'] for y in data['yschema']})==len(data['yschema']),blocks_not_averaged=len(data['schema'])==54 and all(sum(s['receptor_block']==b for s in data['schema'])==27 for b in BLOCKS),
        geometry_na_preserved=all(r['alpha2|gamma2_TYR58__pi_stack__centdist_mean'] is None for r in data['raw']),
        filtered_features_retained_in_raw=all(s['feature_key'] in data['raw'][0] for s in data['audit']),
        mismatch_preserved=all(t['structure_pharmacology_context_match_status']=='context_mismatch' for t in data['targets']),
        source_ids_present=all(t['source_activity_id'] and t['assay_id'] and t['document_id'] and t['source_sha256'] for t in data['targets']),
        missing_y_not_zero=all(r[y['column']] is None for r in data['dataset'] for y in data['yschema'] if not y['availability'].startswith('observed')),
        n_two_no_training=len(data['dataset'])==2 and all(m['state']=='definition_and_features_only; no_training_no_CV_no_importance' for m in data['manifest']),
        baseline_union_bounded=all(0<=v<=1 for r in data['baseline'] for k,v in r.items() if k.endswith('__frequency')))
    if not all(checks.values()):raise ValueError(checks)
    rf.write_json(path/'reports/QC_REPORT.json',dict(status='PASS',checks=checks,tests=test_summary,n_drugs=2,x_features=54,observed_y=4,retained_frequency_features=len(major),retained_geometry_features=len(data['reduced'][0])-1-len(major),training_executed=False))
    report(data,path,config)
    rf.write_json(path/'run_manifest.json',dict(phase='structure_activity_ml_prototype',created_at=datetime.now().astimezone().isoformat(),parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),analysis_contract_sha256=rf.sha((ROOT/'ANALYSIS_CONTRACT.md').read_bytes()),status='completed',artifact_sha256={str(p.relative_to(path)):rf.sha(p.read_bytes()) for p in sorted(path.rglob('*')) if p.is_file()}))
    print(path)


def report(data,path,config):
    kept=[s for s in data['audit'] if s['retained_for_presentation'] and s['feature_type']!='geometry']
    lines=['# Structure–activity ML dataset prototype','', '## Research question',
        'Residue interaction fingerprintに、公開DBのpharmacological activity profileを説明・予測する情報が含まれるか。主解析は **same receptor → different drugs**。本runはX–y形式の構築とQCのみで、学習・CV・feature importanceを実行しない。',
        '', '## Current evidence',
        '薬剤2行、構造Xはα1 27列＋α2 27列（36 frequency/type、18 geometry）。同一構造block内で薬剤差があり、同一記録assay context内でもKiの数値差がある。対応表を構築できたことは予測性能の証拠ではない。',
        '', '| Receptor-fixed pharmacology | Diazepam Ki (nM) | Alprazolam Ki (nM) | Δ(alp−dia), nM | Dia/alp | Assay |', '|---|---:|---:|---:|---:|---|']
    for p in data['pharm']:lines.append(f"| {p['receptor_composition']} | {p['diazepam_value']} | {p['alprazolam_value']} | {p['drug_difference']} | {float(p['diazepam_over_alprazolam_ratio']):.6g} | {p['diazepam_assay_chembl_id']} |")
    lines+=['', 'Kiは小さいほど報告された結合affinityが高い。すべてhuman、relation =、nM、同一document CHEMBL5143601。ただし不確実性推定はなく、有意差・選択性ラベルへ変換しない。activity IDs、assay descriptions、raw source pointersはtables/pharmacology_targets.csvとprovenance/selected_activity_full.jsonlに保持。',
        '', '### Mechanically selected presentation features', '| Feature | Diazepam | Alprazolam | Δ(alp−dia) |','|---|---:|---:|---:|']
    for s in kept:
        k=s['feature_key'];d,a=[r[k] for r in data['raw']];lines.append(f'| {label(k)} | {d:.6f} | {a:.6f} | {a-d:+.6f} |')
    lines+=['', 'α1ではTYR58 hyd、PHE77 πに増加があり、TYR58 π/PやPHE77 hydなどは薬剤間で同頻度でも表示対象となる。α2ではLYS156 hydの増加（2/9→5/9）、TYR58 hydとASN60 hydの増加を保持。これは受容体間比較ではない。残存頻度が高いことと薬剤識別力があることは別である。',
        '', '## Raw / reduced policy',
        'Rawには定義した54列をすべて保持。完全な既存IFP行（geometryのmin/max/median、観測数、pose/source IDsも含む）はprovenance/source_ifp.jsonlに保存。図とrawのgeometry列はcontact-row平均centroid distance / angle / offsetで、観測contact条件付き・pose独立ではない。P/T混合時の単一平均には限界があり、型別頻度と元contact provenanceを併用する。',
        f"縮約は両薬剤frequency <{config['frequency_threshold']}かつ絶対差≤{config['maximum_absolute_difference']}を除外する表示ルール。整数contacting-pose数/pose分母のFractionで境界を厳密評価する。{len(kept)} frequency列と、残存π頻度に紐づく{len(data['reduced'][0])-1-len(kept)} geometry列を保持する。Geometryにfrequency閾値を適用しない。その他{len(data['raw'][0])-len(data['reduced'][0])}列はrawに残る。欠損頻度は除外せず保持。これは統計的有意差でもML用feature selectionでもない。",
        '', '## X–y context and missingness',
        '構造α1はα1β3γ2L full pentamer、α2は9CTJ由来provisional β3/α2/γ2 local C/D/E。薬理はα1/α2β2γ2であり、観測された4つのX–y対応はすべてcontext_mismatch。薬剤間assayのexact_context_matchとは別のQC軸である。未観測yはnot_assessed_missing_activityとし、存在しない測定にexact/partial/mismatchを捏造しない。',
        'wide datasetはraw54 X列＋観測Ki 2列＋未取得endpoint 12列（各blockのKd、other binding、EC50、IC50、Emax、response）を保持。未取得はNA、unit/relation/assayもNAであり、ここでのreceptor/speciesは将来取得の要求context。取得済みdrug片側のみの活動をこの対応表のyへ転用しない。該当15 eligible unpaired recordsも別provenanceに保持する。IC50はassay機序でbindingの場合もあるため、将来もendpoint名だけでfunctional potencyと分類しない。',
        '0は既存監査済みcontactなし。NAは未観測/未定義であり0補完しない。ml_dataset_prototype.csvのqc列・drug ID・ChEMBL IDは予測featureではない。feature_schema/target_schemaとmanifestがX/y列の契約である。CSVではNAを明示、missingness_mask.csvも保持。',
        '', '## Future model comparison manifest',
        'Model0: 各blockの9 archived Vina poseの最小score。Model1: 6対象siteで指定interaction typesのいずれかがあるbinary。Model2: 同じsiteのcontact pose集合の和集合/9（type頻度の和ではない）。Model3: type別frequencyとP/T。Model4: Model3＋geometry。Model5: α1/α2のModel4を順序固定で連結。Model0–4はblock別に定義し、Model5の増分は事前指定single-block baselineと比較する。baseline_features.csvとfeature_set_manifest.csvに実列名を保存。Model1/2は対象6site・選択interaction types内のunionであり全残基網羅ではない。',
        'すべて同じy endpoint、assay/context、評価drug集合・splitで比較する。欠損geometry等の処理は将来training fold内のみで決める。Reduced表示ルールを学習feature選択に流用しない。endpointを平均/統合しない。',
        '', '## Not claimed',
        '- 特定残基がKi差を原因として生じさせる。', '- 現時点の2剤で予測性能がある。', '- すべての薬理差を構造Fingerprintだけで説明できる。', '- 小さいinteraction差すべてに生物学的意味がある。',
        '同方向に動くfeatureがあっても、独立drug n=2、単一pose生成seed、同一化学系列、構造/薬理context不一致により対応情報の有無は判断不能。9 poseを独立drugサンプルと扱わない。',
        '', '## Next step',
        '最優先は同一組成・species・assay機序/endpointで複数追加薬剤のyを確保し、対応する構造Xを作ること。β2対応構造またはβ3比較薬理でcontext差を分離する。十分なdrug数と独立評価集合を確保した後、Model0〜5を同じendpointで比較する。反証可能な問いは「同一contextの未知drugで、type/geometry追加やmulti-receptor化がdocking score/単純contact頻度よりheld-out予測誤差を改善するか」。改善がなければ追加表現の情報価値仮説を支持しない。現時点では拡張の実行可能性は支持、予測上の価値は判断不能。',
        '', '## Provenance / reproduction / QC',
        'ANALYSIS_CONTRACT.mdのsnapshot/hashを保存し、既存runは変更しない。source_inventory.json、selected_activity_full.jsonl、source_ifp.jsonlから元runとChEMBL取得レスポンスへ遡れる。',
        '`MPLBACKEND=Agg .venv/bin/python -m pytest tests --junitxml=<new-path>/tests.xml`',
        '`MPLBACKEND=Agg .venv/bin/python -m src.structure_activity_ml_prototype --tests-xml <new-path>/tests.xml`',
        'テストとQC結果はreports/tests.xml、reports/QC_REPORT.json。図はPNG/PDF、numeric scaleをfrequency・geometry各単位・Kiに分離する。既存runと値を照合し、block非平均・missing保持・endpoint分離・context flags・source IDs・学習未実行を確認。']
    (path/'reports/FINAL_REPORT.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,default=ROOT/'config/fingerprint_reduction_rules.json');parser.add_argument('--tests-xml',type=Path,required=True)
    args=parser.parse_args();run(args.config,args.tests_xml)
