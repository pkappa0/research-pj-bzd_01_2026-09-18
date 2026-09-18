"""Primary axis: fixed receptor -> different drugs. See ANALYSIS_CONTRACT.md."""
from __future__ import annotations
import argparse
from decimal import Decimal
import json
from pathlib import Path
import shutil
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

from . import run_manager
from .bzd_pharmacology_filter import (read_jsonl, write_csv, write_json, encode, sha,
                                      file_inventory, primary_group_key, key_for, tests_summary)

ROOT = Path(__file__).resolve().parents[1]
DRUGS = ('diazepam', 'alprazolam')
BLOCKS = ('alpha1', 'alpha2')
BASE_METRICS = [('frequency', 'interaction_frequency', 'fraction_of_archived_poses')]
PI_METRICS = [('type_P_frequency', 'stacking_type_frequency', 'fraction_of_archived_poses'),
              ('type_T_frequency', 'stacking_type_frequency', 'fraction_of_archived_poses'),
              ('centdist_mean', 'geometry', 'angstrom'), ('angle_mean', 'geometry', 'degree'),
              ('offset_mean', 'geometry', 'angstrom')]


def verify_manifest(path):
    manifest = json.loads((path / 'run_manifest.json').read_text())
    for name, digest in manifest['artifact_sha256'].items():
        if sha((path / name).read_bytes()) != digest:
            raise ValueError('Input artifact changed: ' + name)


def difference(diazepam, alprazolam):
    if diazepam is None or alprazolam is None:
        return None
    a, b = Decimal(str(diazepam)), Decimal(str(alprazolam))
    return str(b - a) if a.is_finite() and b.is_finite() else None


def structural_comparison(features, ids):
    index = {(r['receptor_subtype'], r['drug_id'], r['feature_id']): r for r in features}
    if len(index) != len(features): raise ValueError('Duplicate structural feature')
    rows = []
    for block in BLOCKS:
        bases = [r for r in features if r['receptor_subtype'] == block and r['drug_id'] == 'diazepam']
        for d in bases:
            a = index[(block, 'alprazolam', d['feature_id'])]
            if d['receptor_id'] != a['receptor_id'] or d['actual_residue_number'] != a['actual_residue_number'] or d['residue_chain'] != a['residue_chain']:
                raise ValueError('Within-block receptor/mapping mismatch')
            if d['molecule_chembl_id'] != ids['diazepam'] or a['molecule_chembl_id'] != ids['alprazolam']:
                raise ValueError('Compound ID mismatch in IFP')
            for metric, typ, unit in BASE_METRICS + (PI_METRICS if d['interaction_type'] == 'pi_stack' else []):
                dv, av = d[metric], a[metric]
                rows.append(dict(receptor_block=block, structural_receptor_id=d['receptor_id'],
                    structure_scope=d['structure_scope'], structure_species=d['structure_species'],
                    feature=d['feature_id'] + '__' + metric, site_label=d['site_label'],
                    interaction_type=d['interaction_type'], metric=metric, feature_type=typ, feature_unit=unit,
                    diazepam_value=dv, alprazolam_value=av, drug_difference=difference(dv, av),
                    difference_definition='alprazolam_minus_diazepam',
                    diazepam_molecule_chembl_id=ids['diazepam'], alprazolam_molecule_chembl_id=ids['alprazolam'],
                    diazepam_n_poses=d['n_poses'], alprazolam_n_poses=a['n_poses'],
                    diazepam_n_interacting_poses=d['n_interacting_poses'], alprazolam_n_interacting_poses=a['n_interacting_poses'],
                    diazepam_source_row_ids=d['source_row_ids'], alprazolam_source_row_ids=a['source_row_ids'],
                    diazepam_data_status=d['data_status'], alprazolam_data_status=a['data_status'],
                    diazepam_geometry_n=d.get(metric.replace('_mean', '_n_observed')) if typ == 'geometry' else None,
                    alprazolam_geometry_n=a.get(metric.replace('_mean', '_n_observed')) if typ == 'geometry' else None,
                    residue_chain=d['residue_chain'], actual_residue_number=d['actual_residue_number'],
                    structural_drug_context_match_status='exact_context_match',
                    structural_context_scope='same archived receptor preparation and pose-generation protocol; not independent experiments',
                    missingness='one_or_both_values_missing' if dv is None or av is None else 'both_observed_or_audited_zero'))
    return rows


def bridge_status(pharmacology, structure):
    """Known differences dominate partial context. No beta-family collapsing."""
    reasons = []
    pcomp, scomp = pharmacology.get('receptor_composition'), structure.get('receptor_composition')
    ps, ss = pharmacology.get('species'), structure.get('species')
    if pcomp and scomp and pcomp != scomp: reasons.append('receptor_composition_mismatch')
    if ps and ss and ps != 'unknown' and ss != 'unknown' and ps != ss: reasons.append('species_mismatch')
    if reasons: return 'context_mismatch', reasons
    if not pcomp or not scomp or not ps or not ss or 'unknown' in (ps, ss): reasons.append('missing_identity_metadata')
    if not structure.get('full_construct_verified'): reasons.append('structural_construct_or_splice_form_not_exactly_matched')
    if not pharmacology.get('full_construct_verified'): reasons.append('pharmacology_construct_details_not_verified')
    return ('partial_context_match', reasons) if reasons else ('exact_context_match', ['verified_recorded_identity_agreement'])


def pharmacology_comparison(activities, ids):
    groups = {}
    eligible = []
    for row in activities:
        if row['molecule_chembl_id'] not in set(ids.values()): continue
        if row['drug_id'] not in ids or row['molecule_chembl_id'] != ids[row['drug_id']]: raise ValueError('Compound ID/name assignment changed')
        if row['classification'] != 'PRIMARY' or not row['point_comparison_usable'] or row['review_flags']: continue
        if row['subtype'] not in BLOCKS or row['standard_relation'] != '=': continue
        group = key_for('CG', primary_group_key(row))
        if group != row['comparison_group_id']: raise ValueError('Comparison group changed')
        groups.setdefault(group, []).append(row); eligible.append(row)
    output, used = [], set()
    for gid, rows in sorted(groups.items()):
        drugs = {drug: [r for r in rows if r['molecule_chembl_id'] == ids[drug]] for drug in DRUGS}
        if any(not drugs[d] for d in DRUGS): continue
        if any(len(drugs[d]) != 1 for d in DRUGS): raise ValueError('Multiple observations in a drug/context: no representative or pooling allowed')
        d, a = [drugs[x][0] for x in DRUGS]
        if primary_group_key(d) != primary_group_key(a): raise ValueError('Recorded context differs')
        if d['subtype'] != a['subtype']: raise ValueError('Cross-receptor comparison forbidden in primary analysis')
        dv, av = Decimal(d['standard_value']), Decimal(a['standard_value'])
        ratio = str(dv / av) if d['standard_type'] in {'Ki', 'Kd'} and d['endpoint_class'] == 'binding' and dv > 0 and av > 0 else None
        struct = dict(receptor_composition=d['subtype'] + '_beta3_gamma2', species='human', full_construct_verified=False)
        status, reasons = bridge_status(dict(receptor_composition=d['receptor_composition'], species=d['species'], full_construct_verified=False), struct)
        result = dict(receptor_block=d['subtype'], comparison_id=gid, receptor_composition=d['receptor_composition'], species=d['species'],
            endpoint_class=d['endpoint_class'], endpoint=d['standard_type'], endpoint_detail=d['endpoint_detail'], unit=d['standard_units'], relation='=',
            diazepam_value=d['standard_value'], alprazolam_value=a['standard_value'], drug_difference=str(av-dv),
            difference_definition='alprazolam_minus_diazepam', diazepam_over_alprazolam_ratio=ratio,
            ratio_interpretation='Ki(diazepam)/Ki(alprazolam); smaller Ki is higher reported binding affinity; no potency/efficacy conversion',
            document=d['document_chembl_id'], document_title=d['document_detail_title'], document_doi=d['document_detail_doi'],
            pharmacology_context_match_status='exact_context_match',
            context_match_scope='two drugs: exact recorded assay/document key; unreported conditions not independently verified',
            matched_recorded_context=primary_group_key(d), structure_pharmacology_context_match_status=status,
            structure_pharmacology_context_reasons=reasons,
            structural_receptor_composition=struct['receptor_composition'],
            structural_scope='full_pentamer_alpha1_beta3_gamma2L' if d['subtype']=='alpha1' else 'provisional_local_CDE_beta3_alpha2_gamma2',
            interpretation_limit='No error estimates in selected records; numerical difference is not a significance test')
        for drug, source in [('diazepam',d),('alprazolam',a)]:
            for field in ['activity_id','molecule_chembl_id','assay_chembl_id','target_chembl_id','document_chembl_id','standard_type',
                          'standard_value','standard_relation','standard_units','value','relation','units','receptor_composition',
                          'species','target_organism','assay_detail_assay_organism','assay_detail_description','assay_detail_assay_parameters',
                          'source_file','source_json_pointer','source_sha256']:
                result[drug+'_'+field] = source.get(field)
            used.add(source['activity_id'])
        output.append(result)
    selected = [r for r in eligible if r['activity_id'] in used]
    unpaired = [dict(r, receptor_fixed_gap_reason='no_other_drug_in_same_exact_recorded_context') for r in eligible if r['activity_id'] not in used]
    return output, selected, unpaired


def integrate(structural, pharmacology):
    rows = []
    for p in pharmacology:
        for f in structural:
            if f['receptor_block'] != p['receptor_block']: continue
            rows.append(dict(f, pharmacology_comparison_id=p['comparison_id'],
                pharmacology_receptor_composition=p['receptor_composition'], pharmacology_species=p['species'],
                pharmacology_endpoint=p['endpoint'], pharmacology_endpoint_class=p['endpoint_class'], pharmacology_unit=p['unit'],
                diazepam_pharmacology_value=p['diazepam_value'], alprazolam_pharmacology_value=p['alprazolam_value'],
                pharmacology_drug_difference=p['drug_difference'], document=p['document'],
                diazepam_activity_id=p['diazepam_activity_id'], alprazolam_activity_id=p['alprazolam_activity_id'],
                diazepam_assay_id=p['diazepam_assay_chembl_id'], alprazolam_assay_id=p['alprazolam_assay_chembl_id'],
                pharmacology_context_match_status=p['pharmacology_context_match_status'],
                structure_pharmacology_context_match_status=p['structure_pharmacology_context_match_status'],
                structure_pharmacology_context_reasons=p['structure_pharmacology_context_reasons'],
                interpretation='context-mismatched juxtaposition, not exact matched structure-pharmacology evidence'))
    return rows


def multi_receptor_matrix(structural):
    schema, matrix = [], {d: {'drug_id':d} for d in DRUGS}
    for index, f in enumerate(structural):
        key = f['receptor_block'] + '|' + f['feature']
        schema.append(dict(column_index=index, feature_key=key, receptor_block=f['receptor_block'], feature_type=f['feature_type'], unit=f['feature_unit']))
        for drug in DRUGS: matrix[drug][key] = f[drug+'_value']
    return list(matrix.values()), schema


def figures(directory, structural, pharm):
    directory.mkdir(exist_ok=True); plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    output=[]
    def save(fig,name):
        for ext in ['png','pdf']:
            p=directory/(name+'.'+ext);fig.savefig(p,dpi=170,bbox_inches='tight');output.append(p)
        plt.close(fig)
    for block in BLOCKS:
        rows=[r for r in structural if r['receptor_block']==block and r['feature_type']!='geometry']
        labels=[r['site_label'].replace('gamma2_','γ2 ').replace('alpha_','α ')+' · '+r['interaction_type'].replace('hydrophobic_interaction','hydrophobic').replace('hydrogen_bond','H-bond').replace('halogen_bond','halogen').replace('pi_stack','π')+(' / '+r['metric'][5] if r['metric'].startswith('type_') else '') for r in rows]
        values=np.array([[np.nan if r[d+'_value'] is None else float(r[d+'_value']) for d in DRUGS] for r in rows])
        fig,ax=plt.subplots(figsize=(8,10),layout='constrained');cmap=plt.get_cmap('YlGnBu').with_extremes(bad='#ddd')
        im=ax.imshow(np.ma.masked_invalid(values),vmin=0,vmax=1,cmap=cmap,aspect='auto')
        ax.set_xticks([0,1],['Diazepam','Alprazolam']);ax.set_yticks(range(len(rows)),labels)
        for i in range(len(rows)):
            for j in range(2):
                v=values[i,j];ax.text(j,i,'NA' if np.isnan(v) else f'{v:.2f}',ha='center',va='center',color='white' if v>.6 else '#15252b')
        ax.set_title(f'{block.replace("alpha","α")} fixed structural block\nDifferent drugs on the same archived receptor',pad=14)
        fig.colorbar(im,ax=ax,label='Contacting poses / 9 archived poses',shrink=.8)
        fig.supxlabel(('6HUP · α1β3γ2L full pentamer' if block=='alpha1' else '9CTJ C/D/E · provisional β3/α2/γ2 local construct')+'\nZero = audited no contact; NA = missing evidence. Geometry remains separate.',fontsize=9)
        save(fig,block+'_drug_ifp_heatmap')
    fig,axes=plt.subplots(1,len(pharm),figsize=(12,5),layout='constrained',squeeze=False)
    for ax,p in zip(axes.flat,pharm):
        vals=[float(p[d+'_value']) for d in DRUGS]
        ax.plot([0,1],vals,'-o',color='#315f83',lw=2,markersize=8);ax.set_xticks([0,1],['Diazepam','Alprazolam']);ax.set_xlim(-.4,1.4)
        ax.set_ylim(0,max(vals)*1.4);ax.set_ylabel(p['endpoint']+' ('+p['unit']+')')
        ax.set_title(p['receptor_composition'].replace('alpha','α').replace('beta','β').replace('gamma','γ').replace('_','')+' · human\n'+p['document']+' / '+p['diazepam_assay_chembl_id'])
        for k,v in enumerate(vals):ax.annotate(f'{v:g}',(k,v),xytext=(0,8),textcoords='offset points',ha='center')
        ax.text(.5,.96,'Δ(alp−dia)='+p['drug_difference']+' nM',transform=ax.transAxes,ha='center',va='top')
        ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Receptor-fixed drug pharmacology · same recorded assay per panel',fontsize=15)
    fig.supxlabel('Smaller Ki = higher reported binding affinity; no uncertainty estimates or selectivity categories.\nAssay-to-assay: exact recorded match. Structure-to-pharmacology: β2/β3 context mismatch.',fontsize=10)
    save(fig,'receptor_fixed_pharmacology')
    fig,ax=plt.subplots(figsize=(13,6));ax.set_xlim(0,13);ax.set_ylim(0,6);ax.axis('off')
    ax.text(6.5,5.6,'Primary comparison: fixed receptor → different drugs',ha='center',fontsize=18,weight='bold')
    for y,drug in [(3.8,'Diazepam'),(2.1,'Alprazolam')]:
        ax.text(.15,y+.25,drug,fontsize=13,weight='bold')
        for x,color,text in [(2.1,'#deedf7','α1 fingerprint\nfrequency | type | geometry'),(6.1,'#e4f0df','α2 fingerprint*\nfrequency | type | geometry'),(10.1,'#f3f1eb','…')]:
            width=3.4 if x<10 else 2.4
            ax.add_patch(FancyBboxPatch((x,y-.1),width,1.05,boxstyle='round,pad=0.08',facecolor=color,edgecolor='#668'))
            ax.text(x+width/2,y+.43,text,ha='center',va='center',fontsize=12)
    for x in [3.8,7.8]:
        ax.annotate('',xy=(x,3.66),xytext=(x,3.05),arrowprops={'arrowstyle':'<->','color':'#536574','lw':2})
        ax.text(x+.12,3.33,'alp − dia',fontsize=10,va='center')
    ax.text(6.5,1.25,'drug = [ α1 block | α2 block | … ]  ·  concatenate; do not average blocks',ha='center',fontsize=13)
    ax.text(6.5,.6,'*α2 is a provisional local construct. Mixed units retain their labels and missingness mask.\nBinding / potency / efficacy stay separate. Secondary receptor comparisons are not shown here.',ha='center',fontsize=10)
    save(fig,'multi_receptor_drug_fingerprint_concept')
    return output


def report(structural, pharm, unpaired, qc):
    lines=['# Receptor-fixed drug comparison', '',
        '主解析契約: [ANALYSIS_CONTRACT.md](config/ANALYSIS_CONTRACT.md)。**same receptor → different drugs**。全drug differenceはalprazolam − diazepam。',
        '元integration runは変更せず、α1とα2を独立した構造blockとして再構成した。受容体間差は今回の主解析表・図に含めない。', '',
        '## Pharmacology: receptorを固定した薬剤比較', '',
        '| receptor | study / assay | diazepam Ki | alprazolam Ki | alp − dia | Ki(dia)/Ki(alp) | assay間QC | 構造―薬理QC |',
        '|---|---|---:|---:|---:|---:|---|---|']
    for p in pharm:
        lines.append(f'| {p["receptor_composition"]} | {p["document"]} / {p["diazepam_assay_chembl_id"]} | {p["diazepam_value"]} nM | {p["alprazolam_value"]} nM | {p["drug_difference"]} nM | {float(p["diazepam_over_alprazolam_ratio"]):.6g} | {p["pharmacology_context_match_status"]} | {p["structure_pharmacology_context_match_status"]} |')
    lines += ['', '両比較ともhuman、Ki、nM、relation `=`。同じreceptor composition、assay ID、document ID、記録されたdescription/parameters等の完全一致を要求した。',
        'exact_context_matchは**薬理における2剤間の記録条件**の一致であり、未記載の濃度・反復・実験誤差まで検証した意味ではない。',
        'Kiが小さい方が報告されたbinding affinityは高い。上表では両contextともalprazolamのKiが低い。ただし誤差情報がなく有意性は判断不能。binding差をpotency/efficacy差に読み替えない。',
        'β2薬理とβ3 IFPは明示的なcomposition差があるため、**structure_pharmacology_context_match_status = context_mismatch**。',
        'α1構造は6HUP α1β3γ2L full pentamer。α2は9CTJ C/D/E（β3/α2/γ2）の暫定局所constructで、完全なmatched receptorではない。既知のβ差をpartial matchに弱めない。',
        f'PRIMARYかつpoint比較可能な2剤の薬理は19行。そのうち同条件で2剤が揃う4行を使用し、残る{len(unpaired)}行は `unpaired_eligible_pharmacology.csv` に残した。異なるstudy・assayを組み合わせて対を作らない。',
        '同条件で2剤が揃うfunctional potency / efficacyの比較は今回の適格データでは得られなかった。旧runのdiazepamのみの応答値は薬剤間比較へ流用しない。', '',
        '## α1 block: diazepam vs alprazolam', '',
        '固定contextは6HUP由来の同一prepared receptor。薬剤間でTYR58 hydrophobicは7/9→9/9、PHE77 πは2/9→4/9。PHE77 hydrophobicは両剤7/9で同じ。',
        'TYR58 πは両剤4/9だがcentdist meanは3.9675→4.1775 Å（+0.2100 Å）。geometry差をエネルギー差や有利/不利のscoreとは解釈しない。',
        'PHE77はdiazepamでP=2/9・T=0/9、alprazolamでP=3/9・T=1/9。混合型のangle meanだけに縮約しない。',
        'HIS102 hydrophobicは0/9→2/9、πは両剤1/9。LYS156 H-bondは両剤1/9、hydrophobicは0/9→1/9。SER205 H-bondは両剤2/9。ASN60 halogenは1/9→0/9。', '',
        '## α2 block: diazepam vs alprazolam', '',
        '固定contextは同一の9CTJ-derived C/D/E local construct。TYR58 hydrophobicは3/9→4/9。TYR58/PHE77 πは両剤とも0/9で、π geometryは両剤ともNA。',
        'HIS102共通位置（実残基HIS101）のhydrophobicは0/9→1/9、πは2/9→3/9。LYS156共通位置（実残基LYS155）のhydrophobicは2/9→5/9。',
        'SER205共通位置（実残基SER204）のH-bondは1/9→0/9。ASN60 hydrophobicは3/9→4/9。H-bond・hydrophobic・halogenは別featureとして保持した。', '',
        '## 全featureのdrug difference', '',
        '各blockは12 residue×interaction頻度、6 P/T頻度、9 geometry mean = 27 feature。2 blocksで54 feature。',
        'frequencyの単位はfraction、geometryはÅ/degree。geometryは観測raw contactに条件づけたmeanで、pose重複のある独立でない観測。欠損は差もNA。', '',
        '| block | feature | diazepam | alprazolam | alp − dia | unit |', '|---|---|---:|---:|---:|---|']
    fmt=lambda x:'NA' if x is None else f'{float(x):.6g}'
    for r in structural:
        lines.append('| '+' | '.join([r['receptor_block'],r['feature'],fmt(r['diazepam_value']),fmt(r['alprazolam_value']),fmt(r['drug_difference']),r['feature_unit']])+' |')
    lines += ['', '## Supported observation', '',
        '- 同じ受容体block内でも薬剤ごとに残基/type/geometry profileは異なる。異ならないfeatureも保存し、都合のよい差だけを選ばない。',
        '- 各β2薬理contextではalprazolamのKiがdiazepamより低い。α1β2γ2の差−13.2 nMとα2β2γ2の差−19.4 nMは別の薬剤間比較として保持する。',
        '- 同じ頻度でもgeometryやP/T構成が異なるため、frequencyだけへの縮約は構造情報を失う。geometryの薬理説明への増分寄与はまだ検証していない。', '',
        '## Suggestive pattern', '',
        '- α1 blockではalprazolamのTYR58 hydrophobic、PHE77 π等が多く、別のβ2薬理contextで観測される低いKiと並置できる。ただしcompositionが一致しないため、対応を検証したとは言えない。',
        '- α2 blockではalprazolamのHIS π、LYS hydrophobic等が多いが、SER H-bondは少ない。単一の「接触の多さ」scoreにはまとめない。β2薬理との関連は未検証の仮説に留まる。', '',
        '## Not supported', '',
        '- featureがaffinity/efficacyを決定すること、selectivityを予測できること。2剤・context mismatchを含む並置から因果や予測性能を主張できない。',
        '- 薬理の薬剤間ratioを受容体間selectivityと呼ぶこと、bindingをpotency/efficacyへ読み替えること。',
        '- poseを独立実験として扱うこと、統合表の54行を54組の独立薬理観測として扱うこと、相関や回帰をここで計算すること。', '',
        '## Multi-receptor representationと次段階', '',
        '`drug = [α1 fingerprint | α2 fingerprint | ...]`。`multi_receptor_drug_fingerprint.csv` は2剤×54featureを固定順で連結し、`feature_schema.csv` にblock/type/unit/order、`multi_receptor_missingness.csv` に欠損maskを保存。',
        'geometryとfrequencyは異なるunitのまま保存し、一つの色スケールで表示・平均・自動正規化しない。連結は表現設計でありML学習ではない。',
        '次に検証すべき仮説は、同一の検証済み受容体・同じ調製条件における**薬剤間feature差**が独立seedでも再現するか。',
        '例：α1固定でalprazolam−diazepamのPHE77 π頻度差が正という仮説。条件/pose採択を事前固定した5 seedsの4以上で正、かつseed中央値が正を暫定基準とし、満たさなければこのfeature差の頑健性を支持しない。',
        'その後、同じβ背景・同じ受容体compositionの両剤Kiと誤差を取得して検証する。α2も別block内のdrug differenceとして検証し、主解析を受容体間差へ戻さない。', '',
        '## 図', '',
        '![α1 fixed block](figures/alpha1_drug_ifp_heatmap.png)',
        '![α2 fixed block](figures/alpha2_drug_ifp_heatmap.png)',
        '![Fixed receptor pharmacology](figures/receptor_fixed_pharmacology.png)',
        '![Multi-receptor concept](figures/multi_receptor_drug_fingerprint_concept.png)', '',
        '## QC・再現', '',
        f'QC {qc["status"]}: {len(qc["checks"])} checks PASS。tests `{encode(qc["tests"])}`。',
        '元integration/pharmacology runの全artifact hashを検証。元IFP数値、薬理4行全field、ID・単位・relationを保持。元データと過去runは変更しない。',
        'contract、解析コード、tests、入力snapshotのSHA-256をmanifestへ記録。Git commitは実行前HEAD、実行コードはsnapshot hashで特定。LATEST_RUN.txtは不変。',
        '再実行: `MPLBACKEND=Agg python -m pytest -q --junitxml=/tmp/receptor-fixed-tests.xml` 後、',
        '`python -m src.receptor_fixed_drug_comparison --test-results /tmp/receptor-fixed-tests.xml`。', '']
    return '\n'.join(lines)


def prepare(config):
    ifp=ROOT/config['ifp_run'];pharm=ROOT/config['pharmacology_run']
    verify_manifest(ifp);verify_manifest(pharm)
    features=read_jsonl(ifp/'processed/diazepam_alprazolam_ifp.jsonl')
    activities=read_jsonl(pharm/'processed/bzd_activity_classification.jsonl')
    structural=structural_comparison(features,config['compound_ids'])
    pharmacology,selected,unpaired=pharmacology_comparison(activities,config['compound_ids'])
    integrated=integrate(structural,pharmacology);matrix,schema=multi_receptor_matrix(structural)
    checks={
        'primary_axis_is_within_receptor_drug_difference':all(r['difference_definition']=='alprazolam_minus_diazepam' for r in structural+pharmacology),
        'two_blocks_54_distinct_features':len(structural)==54 and len({(r['receptor_block'],r['feature']) for r in structural})==54,
        'all_drug_differences_correct':all(r['drug_difference']==difference(r['diazepam_value'],r['alprazolam_value']) for r in structural+pharmacology),
        'geometry_missing_not_zero':all(r['drug_difference'] is None for r in structural if r['diazepam_value'] is None or r['alprazolam_value'] is None),
        'two_fixed_assay_comparisons_four_source_ids':len(pharmacology)==2 and len({r['activity_id'] for r in selected})==4,
        'both_drugs_same_assay_document_composition_species_endpoint_unit_relation':all(all(p['diazepam_'+k]==p['alprazolam_'+k] for k in ['assay_chembl_id','document_chembl_id','receptor_composition','species','standard_type','standard_units','standard_relation']) for p in pharmacology),
        'beta2_beta3_never_matched':all(p['structure_pharmacology_context_match_status']=='context_mismatch' for p in pharmacology),
        'integrated_receptor_blocks_not_crossed':all(r['receptor_block'] in r['pharmacology_receptor_composition'] for r in integrated),
        'no_activity_value_modification':all(r==next(a for a in activities if a['activity_id']==r['activity_id']) for r in selected),
        'all_19_eligible_rows_retained_once':len(selected)+len(unpaired)==19 and len({r['activity_id'] for r in selected+unpaired})==19,
        'multi_receptor_blocks_preserved':len(matrix)==2 and len(schema)==54 and len({r['feature_key'] for r in schema})==54,
    }
    if not all(checks.values()):raise ValueError('QC failed: '+str(checks))
    return structural,pharmacology,selected,unpaired,integrated,matrix,schema,checks


def tracked_inventory():
    names=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
    return file_inventory(ROOT,[ROOT/n for n in names if n])


def run(config_path,test_results=None):
    config=json.loads(config_path.read_text());before=tracked_inventory()
    contract=ROOT/config['analysis_contract']
    if 'same receptor → different drugs' not in contract.read_text():raise ValueError('Missing primary analysis contract')
    structural,pharm,selected,unpaired,integrated,matrix,schema,checks=prepare(config)
    source_names=[config['ifp_run']+'/processed/diazepam_alprazolam_ifp.jsonl',config['ifp_run']+'/pi_stacking_geometry_contacts.csv',
                  config['ifp_run']+'/run_manifest.json',config['pharmacology_run']+'/processed/bzd_activity_classification.jsonl',config['pharmacology_run']+'/run_manifest.json']
    code_names=['src/receptor_fixed_drug_comparison.py','src/bzd_pharmacology_filter.py','src/acquire_chembl.py','src/run_manager.py','tests/test_receptor_fixed_drug_comparison.py','pytest.ini']
    hashes=file_inventory(ROOT,[ROOT/n for n in source_names]);test_info=tests_summary(test_results)
    out,manifest=run_manager.create_run(short_phase_name='receptor_fixed_drug_comparison',analysis_phase='Same receptor, different drugs',
        analysis_purpose='Describe within-receptor drug feature contrasts and context-qualified pharmacology',parent_run=Path(config['ifp_run']).name,
        drugs=list(DRUGS),receptors={'alpha1':'6HUP alpha1 beta3 gamma2L','alpha2':'9CTJ provisional CDE local construct'},docking_parameters={'new_docking':False},
        input_files=source_names,notes=['Primary axis: same receptor -> different drugs','All differences alprazolam minus diazepam','No changes to prior runs or LATEST_RUN.txt'],code_paths=[ROOT/n for n in code_names])
    try:
        for n in source_names:
            dest=out/'raw/lineage'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,dest)
        for n in code_names:
            dest=out/'config/code_snapshot'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,dest)
        shutil.copyfile(contract,out/'config/ANALYSIS_CONTRACT.md');shutil.copyfile(config_path,out/'config/receptor_fixed_drug_comparison.json')
        tables={'receptor_fixed_drug_comparison':structural,'receptor_fixed_pharmacology':pharm,'structure_pharmacology_drug_comparison':integrated,
                'unpaired_eligible_pharmacology':unpaired,'multi_receptor_drug_fingerprint':matrix,'feature_schema':schema,
                'multi_receptor_missingness':[{k:(v is None) if k!='drug_id' else v for k,v in r.items()} for r in matrix]}
        for name,records in tables.items():write_csv(out/(name+'.csv'),records)
        for name,records in [('selected_pharmacology',selected),('receptor_fixed_drug_comparison',structural),('structure_pharmacology_drug_comparison',integrated)]:
            (out/'processed'/(name+'.jsonl')).write_text(''.join(encode(r)+'\n' for r in records))
        if test_results:shutil.copyfile(test_results,out/'logs/pytest.xml')
        images=figures(out/'figures',structural,pharm)
        checks['four_figures_png_pdf_nonempty']=len(images)==8 and all(p.stat().st_size>1000 for p in images)
        checks['source_snapshots_identical']=all(sha((out/'raw/lineage'/n).read_bytes())==h for n,h in hashes.items())
        checks['contract_snapshot_identical']=sha(contract.read_bytes())==sha((out/'config/ANALYSIS_CONTRACT.md').read_bytes())
        checks['preexisting_tracked_files_unchanged']=tracked_inventory()==before
        checks['null_roundtrip']=read_jsonl(out/'processed/receptor_fixed_drug_comparison.jsonl')==structural
        if not all(checks.values()):raise ValueError('Artifact QC failed')
        qc=dict(status='PASS',checks=checks,tests=test_info,source_input_sha256=hashes,preexisting_tracked_file_count=len(before),
                structural_feature_count=len(structural),pharmacology_comparison_count=len(pharm),pharmacology_source_activity_count=len(selected),
                structure_pharmacology_exact_match_count=0,analysis_contract_sha256=sha(contract.read_bytes()))
        write_json(out/'QC_REPORT.json',qc);write_json(out/'reports/preexisting_files_sha256.json',before)
        (out/'FINAL_REPORT.md').write_text(report(structural,pharm,unpaired,qc))
        manifest.update(completed_at=run_manager.iso_now(),qc_status='PASS_WITH_CONTEXT_MISMATCH',analysis_contract_sha256=qc['analysis_contract_sha256'],source_input_sha256=hashes,
            artifact_sha256=file_inventory(out,[p for p in out.rglob('*') if p.is_file() and p!=out/'run_manifest.json']))
        manifest['output_files']=list(manifest['artifact_sha256']);run_manager.write_manifest(out,manifest)
    except Exception as exc:
        manifest.update(completed_at=run_manager.iso_now(),qc_status='FAILED',error=str(exc));run_manager.write_manifest(out,manifest);raise
    return out,qc


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=ROOT/'config/receptor_fixed_drug_comparison.json');p.add_argument('--test-results',type=Path)
    a=p.parse_args();out,qc=run(a.config,a.test_results);print(encode({'run_directory':str(out.relative_to(ROOT)),'qc':qc['status']}))


if __name__=='__main__':main()
