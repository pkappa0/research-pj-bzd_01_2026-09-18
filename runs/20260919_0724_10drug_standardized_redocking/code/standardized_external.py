"""External annotation only after verifying the structural freeze."""
from .standardized_redocking import *
from collections import defaultdict,Counter
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    run=runpath();freeze=json.loads((run/'STRUCTURAL_RESULTS_FROZEN.json').read_text());assert all(sha(run/p)==h for p,h in freeze['artifacts'].items())
    assert sha(run/'config/unsupervised_analysis_config.json')==freeze['config_sha256']
    started=datetime.now().astimezone().isoformat();assert started>freeze['frozen_at']
    source=ROOT/'runs/20260918_1518_bzd_pharmacology_filter/processed/bzd_activity_classification.jsonl'
    activities=[json.loads(l) for l in source.read_text().splitlines()]
    eligible=[r for r in activities if r['drug_id'] in DRUGS and r['classification']=='PRIMARY' and r['point_comparison_usable'] and not r['review_flags'] and r['standard_relation']=='=']
    groups=defaultdict(list)
    for r in eligible:groups[r['comparison_group_id']].append(r)
    matched={g:rr for g,rr in groups.items() if len({r['drug_id'] for r in rr})>=2 and len({r['drug_id'] for r in rr})==len(rr)}
    annotation=[];pairs=[];correspond=[]
    distances={rep:{r['drug_id']:r for r in read(run/'tables/similarity_euclidean.csv') if r['representation']==rep} for rep in ['B_alpha1','B_alpha2','E']}
    clusters={rep:{r['drug_id']:r['cluster'] for r in read(run/'tables/cluster_assignments.csv') if r['representation']==rep} for rep in distances}
    for gid,rows in matched.items():
        for r in rows:
            annotation.append(dict(drug_id=r['drug_id'],molecule_chembl_id=r['molecule_chembl_id'],receptor=r['subtype'],endpoint_class=r['endpoint_class'],endpoint=r['standard_type'],value=r['standard_value'],unit=r['standard_units'],relation=r['standard_relation'],receptor_composition=r['receptor_composition'],species=r['species'],assay_id=r['assay_chembl_id'],target_id=r['target_chembl_id'],document_id=r['document_chembl_id'],activity_id=r['activity_id'],comparison_group_id=gid,assay_description=r['assay_detail_description'],assay_parameters=r['assay_detail_assay_parameters'],source_file=r['source_file'],source_sha256=r['source_sha256'],source_json_pointer=r['source_json_pointer'],source_run=str(source.relative_to(ROOT)),pharmacology_match='exact_recorded_context',structure_pharmacology_match='context_mismatch' if 'beta2' in r['receptor_composition'] else 'partial_context_match',availability='observed_matched'))
        for a,b in itertools.combinations(rows,2):
            d,e=a['drug_id'],b['drug_id'];rep='B_'+a['subtype'];pairs.append(dict(kind='pharmacology',drug_a=d,drug_b=e,context_id=gid,receptor=a['subtype'],endpoint=a['standard_type'],unit=a['standard_units'],value_a=a['standard_value'],value_b=b['standard_value'],difference_b_minus_a=float(b['standard_value'])-float(a['standard_value']),structural_distance=float(distances[rep][d][e]),same_primary_cluster=clusters[rep][d]==clusters[rep][e],same_E_cluster=clusters['E'][d]==clusters['E'][e],activity_id_a=a['activity_id'],activity_id_b=b['activity_id'],assay_id=a['assay_chembl_id'],document_id=a['document_chembl_id'],context_status='recorded_same_assay; structure_beta3_vs_pharmacology_beta2_mismatch',interpretation='descriptive pair only; no pooled endpoint or association test'))
    for d in DRUGS:
        for b in BLOCKS:
            for family in ['binding','functional_potency','functional_efficacy']:
                if not any(r['drug_id']==d and r['receptor']==b and r['endpoint_class']==family for r in annotation):annotation.append(dict(drug_id=d,receptor=b,endpoint_class=family,value=None,endpoint=None,availability='no_matched_recorded_context_in_existing_filter',structure_pharmacology_match='not_assessed_missing_activity'))
    table(run/'tables/external_pharmacology_annotation.csv',annotation)
    selected_ids={r['activity_id'] for rr in matched.values() for r in rr}
    (run/'logs/external_selected_activity_full.jsonl').write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in eligible if r['activity_id'] in selected_ids))
    (run/'logs/external_unpaired_activity_full.jsonl').write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in eligible if r['activity_id'] not in selected_ids))
    phen_source=ROOT/'runs/20260917_1706_structure_mouse_bridge/data/processed/mouse_phenotype_raw.csv'
    phen=[dict(r,source_row_id=i,source_file=str(phen_source.relative_to(ROOT)),source_sha256=sha(phen_source)) for i,r in enumerate(read(phen_source)) if r['species']=='mouse' and r['evidence_tier']=='P1']
    pg=defaultdict(list);condition=['study_id','species','strain','assay','route','dose_mg_per_kg','time_min','endpoint','raw_unit']
    for r in phen:
        if all(r[k] and r[k]!='not_reported_in_extracted_table' for k in ['dose_mg_per_kg','time_min','raw_value']):pg[tuple(r[k] for k in condition)].append(r)
    matchedphen={k:v for k,v in pg.items() if len({r['drug_id'] for r in v})>=2};matchedrows={r['source_row_id'] for rr in matchedphen.values() for r in rr}
    for r in phen:r['comparison_status']='same_recorded_dose_time_context; strain_unreported' if r['source_row_id'] in matchedrows else 'context_separate_not_pooled'
    for d in DRUGS:
        if not any(r['drug_id']==d for r in phen):phen.append(dict(drug_id=d,raw_value=None,comparison_status='no_existing_P1_mouse_evidence',availability='NA'))
    table(run/'tables/external_phenotype_annotation.csv',phen)
    for key,rows in matchedphen.items():
        for a,b in itertools.combinations(rows,2):
            d,e=a['drug_id'],b['drug_id'];pairs.append(dict(kind='mouse_phenotype',drug_a=d,drug_b=e,context_id='|'.join(key),receptor='E',endpoint=a['endpoint'],unit=a['raw_unit'],value_a=a['raw_value'],value_b=b['raw_value'],difference_b_minus_a=float(b['raw_value'])-float(a['raw_value']),structural_distance=float(distances['E'][d][e]),same_primary_cluster=clusters['E'][d]==clusters['E'][e],study_id=a['study_id'],dose_mg_per_kg=a['dose_mg_per_kg'],time_min=a['time_min'],route=a['route'],source_row_a=a['source_row_id'],source_row_b=b['source_row_id'],context_status='same_recorded_conditions; strain_missing; repeated_conditions_not_independent_drug_pairs',interpretation='counts and conditions retained; not normalized to a potency label'))
    table(run/'tables/external_matched_comparisons.csv',pairs)
    from rdkit import Chem
    compounds={r['drug_id']:r for r in read(ROOT/'runs/20260918_1117_chembl_acquisition/processed/compounds.csv')}
    for d in DRUGS:
        m=Chem.MolFromSmiles(compounds[d]['canonical_smiles']);benz=any(len(r)==7 and sum(m.GetAtomWithIdx(i).GetSymbol()=='N' for i in r)==2 for r in m.GetRingInfo().AtomRings())
        nearest=min((e for e in DRUGS if e!=d),key=lambda e:(float(distances['E'][d][e]),e))
        data=[r for r in annotation if r['drug_id']==d and r['availability']=='observed_matched'];ph=[r for r in phen if r['drug_id']==d and r.get('evidence_tier')=='P1']
        correspond.append(dict(drug_id=d,alpha1_cluster=clusters['B_alpha1'][d],alpha2_cluster=clusters['B_alpha2'][d],multireceptor_cluster=clusters['E'][d],nearest_E_drug=nearest,nearest_E_distance=float(distances['E'][d][nearest]),postfreeze_chemical_class='BZD-like_7ring_2N' if benz else 'non-BZD_ring_topology',class_source='post-freeze RDKit ring topology on archived ChEMBL SMILES; descriptive chemical annotation only',matched_pharmacology_records=len(data),matched_pharmacology_profile=[{k:r[k] for k in ['receptor','endpoint','value','unit','assay_id','document_id']} for r in data] or None,mouse_P1_records=len(ph),mouse_studies=sorted({r['study_id'] for r in ph}) or None,functional_profile=None,assessment='external_coverage_insufficient_for_general_correspondence'))
    table(run/'tables/structure_external_correspondence.csv',correspond)
    fig,ax=plt.subplots(figsize=(19,8));ax.axis('off')
    cells=[]
    for r in correspond:
        d=r['drug_id'];p=r['matched_pharmacology_profile'];lines=[]
        if p:
            for t in p:lines.append(t['receptor']+' '+t['endpoint']+' '+str(float(t['value']))+' '+t['unit']+' ['+t['assay_id'].replace('CHEMBL','')+']')
        cells.append([d,r['alpha1_cluster'],r['alpha2_cluster'],r['multireceptor_cluster'],r['postfreeze_chemical_class'].replace('_',' '),'\n'.join(lines) if lines else 'NA', '\n'.join(r['mouse_studies']) if r['mouse_studies'] else 'NA'])
    t=ax.table(cellText=cells,colLabels=['Drug','α1 C','α2 C','E C','Post-freeze class','Matched binding (assay IDs)','Mouse P1 studies'],cellLoc='left',loc='center',colWidths=[.12,.04,.04,.04,.19,.36,.16]);t.auto_set_font_size(False);t.set_fontsize(8)
    for (i,j),cell in t.get_celld().items():cell.set_height(.064 if i==0 else .073 if i!=1 else .19)
    ax.set_title('Frozen structural clusters → external annotation (no refitting)\nKi comparisons are assay-specific; β2 pharmacology / β3 structure mismatch; matched functional activity NA',fontsize=14)
    fig.text(.05,.02,'Mouse observations retain dose, time, assay and study; study availability is not a phenotype score. Cluster numbers are arbitrary within each representation.',fontsize=10)
    for ext in ['png','pdf']:fig.savefig(run/'figures'/('structure_pharmacology_phenotype_overview.'+ext),dpi=150,bbox_inches='tight')
    plt.close(fig)
    save(run/'reports/EXTERNAL_ANNOTATION_PROVENANCE.json',dict(started_at=started,structural_freeze_at=freeze['frozen_at'],freeze_sha256=sha(run/'STRUCTURAL_RESULTS_FROZEN.json'),sources={str(source.relative_to(ROOT)):sha(source),str(phen_source.relative_to(ROOT)):sha(phen_source)},matched_pharmacology_contexts=len(matched),matched_activity_count=len(selected_ids),mouse_P1_rows=sum(r.get('evidence_tier')=='P1' for r in phen),matched_mouse_conditions=len(matchedphen),independent_mouse_drug_pairs=len({tuple(sorted((r['drug_a'],r['drug_b']))) for r in pairs if r['kind']=='mouse_phenotype'})))
    assert all(sha(run/p)==h for p,h in freeze['artifacts'].items())
    print('ANNOTATION COMPLETE',len(matched),'pharmacology contexts',len(matchedphen),'mouse condition pairs',flush=True)

if __name__=='__main__':main()
