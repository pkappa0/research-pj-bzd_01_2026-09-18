import json
from copy import deepcopy
from decimal import Decimal
import pytest
from src import receptor_fixed_drug_comparison as rf

@pytest.fixture(scope='module')
def data():
    config=json.loads((rf.ROOT/'config/receptor_fixed_drug_comparison.json').read_text())
    return config,rf.prepare(config)

def test_drug_axis_not_receptor_axis(data):
    rows=data[1][0]
    assert len(rows)==54
    r=next(r for r in rows if r['receptor_block']=='alpha1' and r['feature']=='gamma2_TYR58__hydrophobic_interaction__frequency')
    assert r['diazepam_value']==7/9 and r['alprazolam_value']==1
    assert float(r['drug_difference'])==pytest.approx(2/9)
    r=next(r for r in rows if r['receptor_block']=='alpha2' and r['feature']=='alpha_LYS156__hydrophobic_interaction__frequency')
    assert float(r['drug_difference'])==pytest.approx(3/9)

def test_geometry_differences_and_missing(data):
    rows=data[1][0]
    r=next(r for r in rows if r['receptor_block']=='alpha1' and r['feature']=='gamma2_TYR58__pi_stack__centdist_mean')
    assert float(r['drug_difference'])==pytest.approx(.21)
    r=next(r for r in rows if r['receptor_block']=='alpha2' and r['feature']=='gamma2_TYR58__pi_stack__centdist_mean')
    assert r['diazepam_value'] is None and r['alprazolam_value'] is None and r['drug_difference'] is None

def test_drug_pharmacology_exact_recorded_match_but_bridge_mismatch(data):
    p=data[1][1]
    assert len(p)==2 and {r['drug_difference'] for r in p}=={'-13.2','-19.4'}
    assert all(r['pharmacology_context_match_status']=='exact_context_match' for r in p)
    assert all(r['structure_pharmacology_context_match_status']=='context_mismatch' for r in p)
    assert all(r['diazepam_assay_chembl_id']==r['alprazolam_assay_chembl_id'] for r in p)

@pytest.mark.parametrize('change,expected',[({},'exact_context_match'),({'full_construct_verified':False},'partial_context_match'),({'receptor_composition':'alpha1_beta3_gamma2'},'context_mismatch'),({'species':'rat'},'context_mismatch'),({'species':None},'partial_context_match')])
def test_context_status_scopes(change,expected):
    p=dict(receptor_composition='alpha1_beta2_gamma2',species='human',full_construct_verified=True)
    s=dict(p,**{}) ; s.update(change)
    assert rf.bridge_status(p,s)[0]==expected

def test_mismatch_not_weakened_by_local_construct():
    p=dict(receptor_composition='alpha2_beta2_gamma2',species='human')
    s=dict(receptor_composition='alpha2_beta3_gamma2',species='human',full_construct_verified=False)
    assert rf.bridge_status(p,s)[0]=='context_mismatch'

def test_no_cross_study_pair_or_duplicate_representative(data):
    config=data[0];rows=rf.read_jsonl(rf.ROOT/config['pharmacology_run']/'processed/bzd_activity_classification.jsonl')
    rows.append(deepcopy(next(r for r in rows if r['activity_id']==24843226)))
    with pytest.raises(ValueError,match='Multiple observations'):rf.pharmacology_comparison(rows,config['compound_ids'])

def test_matrix_retains_blocks_units_and_missing(data):
    matrix,schema=data[1][5:7]
    assert len(matrix)==2 and len(schema)==54
    assert sum(r['receptor_block']=='alpha1' for r in schema)==27
    assert sum(r['receptor_block']=='alpha2' for r in schema)==27
    assert matrix[0]['alpha2|gamma2_TYR58__pi_stack__centdist_mean'] is None
    assert {r['unit'] for r in schema}=={'fraction_of_archived_poses','angstrom','degree'}

def test_integrated_rows_are_context_qualified_not_replicates(data):
    rows=data[1][4]
    assert len(rows)==54 and len({r['pharmacology_comparison_id'] for r in rows})==2
    assert all(r['receptor_block'] in r['pharmacology_receptor_composition'] for r in rows)
    assert all(r['structure_pharmacology_context_match_status']=='context_mismatch' for r in rows)

def test_preserve_unpaired_and_source_fields(data):
    selected,unpaired=data[1][2:4]
    assert len(selected)==4 and len(unpaired)==15
    assert len({r['activity_id'] for r in selected+unpaired})==19
    assert all(r['standard_relation']=='=' for r in selected)

def test_contract_and_report_axis(data):
    text=(rf.ROOT/'ANALYSIS_CONTRACT.md').read_text()
    assert 'same receptor → different drugs' in text
    assert 'same drug → different receptors' in text
    assert 'drug = [α1 fingerprint | α2 fingerprint | ...]' in text
    assert 'binding / potency / efficacy' in text

def test_missing_never_becomes_zero():
    assert rf.difference(None,None) is None
    assert rf.difference(None,1) is None
    assert rf.difference('14.0','0.8')=='-13.2'
