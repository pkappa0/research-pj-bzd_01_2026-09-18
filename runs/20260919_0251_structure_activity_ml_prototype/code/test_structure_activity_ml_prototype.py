import json
from fractions import Fraction as F
import pytest
from src import structure_activity_ml_prototype as ml

@pytest.fixture(scope='module')
def config():return json.loads((ml.ROOT/'config/fingerprint_reduction_rules.json').read_text())
@pytest.fixture(scope='module')
def data(config):return ml.prepare(config)

@pytest.mark.parametrize('d,a,keep',[(F(0),F(2,9),False),(F(1,9),F(3,9),False),(F(0),F(3,9),True),(F('0.4'),F('0.4'),True),(None,F(0),True)])
def test_threshold_boundary(d,a,keep,config):assert ml.display_keep(d,a,config)==keep

def test_frequency_validation():
    with pytest.raises(ValueError):ml.frequency_fraction(dict(n_interacting_poses=2,n_poses=9,frequency=.5),'frequency')

def test_raw_source_exact(data):
    for row in data['structural']:
        for drug in ml.DRUGS:
            s=next(s for s in data['inputs'] if s['drug_id']==drug and s['receptor_subtype']==row['receptor_block'] and s['feature_id']==row['feature'].rsplit('__',1)[0])
            assert row[drug+'_value']==s[row['metric']]
    assert len(data['schema'])==54

def test_reduction_and_geometry(data):
    retained=[s for s in data['audit'] if s['retained_for_presentation']]
    assert len(retained)==15
    assert sum(s['feature_type']=='geometry' for s in retained)==6
    assert len(data['raw'][0])==55 and len(data['reduced'][0])==16
    assert all(s['feature_key'] in data['raw'][0] for s in data['audit'])

def test_observed_y_and_context(data):
    expected={24843214:'14.0',24843226:'0.8',24843235:'20.0',24843247:'0.6'}
    assert {r['source_activity_id']:r['value'] for r in data['targets']}==expected
    assert all(r['pharmacology_context_match_status']=='exact_context_match' and r['structure_pharmacology_context_match_status']=='context_mismatch' for r in data['targets'])
    assert {p['drug_difference'] for p in data['pharm']}=={'-13.2','-19.4'}
    assert len(data['unpaired'])==15

def test_missing_y_and_geometry(data):
    assert len(data['yschema'])==14
    for r in data['dataset']:
        assert r['alpha2|gamma2_TYR58__pi_stack__centdist_mean'] is None
        assert r['alpha2|gamma2_TYR58__pi_stack__frequency']==0
        assert sum(r[s['column']] is None for s in data['yschema'])==12
    assert all(m['context_match_status']=='not_assessed_missing_activity' for m in data['mappings'] if m['value'] is None)

def test_manifest_no_fitting_and_real_columns(data):
    assert len(data['manifest'])==11 and len(data['dataset'])==2
    for m in data['manifest']:
        table=data['raw'] if m['feature_table']=='raw_fingerprint.csv' else data['baseline']
        assert all(k in table[0] for k in m['columns'])
        assert 'no_training_no_CV_no_importance' in m['state']
    assert next(m for m in data['manifest'] if m['model']==5)['feature_count']==54

def test_type_collapsed_frequency_is_union(data):
    for drug,r in zip(ml.DRUGS,data['baseline']):
        for block in ml.BLOCKS:
            for site in {'gamma2_TYR58','gamma2_PHE77','alpha_HIS102','alpha_LYS156','alpha_SER205','gamma2_ASN60'}:
                source=[s for s in data['inputs'] if s['drug_id']==drug and s['receptor_subtype']==block and s['site_label']==site]
                union=set(p for s in source for p in s['pose_ids'])
                assert r[f'{block}|{site}__any_selected_contact__frequency']==len(union)/9

def test_csv_explicit_na_roundtrip(tmp_path):
    p=tmp_path/'x.csv';ml.csv_write(p,[{'drug':'a','zero':0,'missing':None}])
    assert p.read_text().splitlines()[1]=='a,0,NA'
