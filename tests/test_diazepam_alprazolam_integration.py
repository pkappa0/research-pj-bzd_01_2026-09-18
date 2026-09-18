"""Numerical comparison guards, raw IFP provenance and immutable integration run tests."""
from copy import deepcopy
from decimal import Decimal
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from src import diazepam_alprazolam_integration as it


@pytest.fixture(scope='module')
def real():
    config = json.loads((it.ROOT / 'config/diazepam_alprazolam_integration.json').read_text())
    return config, it.prepare(config)


def numeric(value, **kwargs):
    return dict(standard_type='Ki', standard_units='nM', endpoint_class='binding', standard_relation='=', standard_value=value, **kwargs)


def test_exact_decimal_difference_and_ki_ratio():
    r = it.arithmetic(numeric('0.8'), numeric('0.6'))
    assert r['alpha2_minus_alpha1'] == '-0.2'
    assert Decimal(r['alpha1_over_alpha2_ratio']) == Decimal('0.8') / Decimal('0.6')
    assert r['ratio_eligible']


@pytest.mark.parametrize('field,value', [('standard_relation','>'),('standard_relation','<'),('standard_units','uM'),
    ('standard_units',None),('standard_type','Kd'),('endpoint_class','functional_potency'),('standard_value',None),
    ('standard_value','NaN'),('standard_value','Infinity')])
def test_incompatible_or_censored_data_never_get_ratio(field, value):
    a, b = numeric('14'), numeric('20'); b[field] = value
    r = it.arithmetic(a, b)
    assert r['alpha1_over_alpha2_ratio'] is None and not r['ratio_eligible']
    assert r['alpha2_minus_alpha1'] is None


@pytest.mark.parametrize('value', ['0', '-1'])
def test_nonpositive_binding_ratio_withheld(value):
    r = it.arithmetic(numeric('14'), numeric(value))
    assert not r['ratio_eligible']


def test_percent_change_is_not_automatically_ratio_scale():
    a, b = numeric('156.0'), numeric('89.0')
    for r in (a, b): r.update(standard_type='Efficacy', standard_units='%', endpoint_class='functional_efficacy')
    r = it.arithmetic(a, b)
    assert r['alpha2_minus_alpha1'] == '-67.0'
    assert r['difference_unit'] == 'percentage_points'
    assert r['alpha1_over_alpha2_ratio'] is None


def test_real_pairs_not_averaged(real):
    config, data = real
    pairs, activities = data[:2]
    assert len(pairs) == 4 and len(activities) == 8
    ki = {(p['drug_id'], p['document_chembl_id']): p for p in pairs if p['standard_type'] == 'Ki'}
    assert ki[('diazepam', 'CHEMBL5143601')]['alpha2_minus_alpha1'] == '6.0'
    assert ki[('diazepam', 'CHEMBL6078658')]['alpha2_minus_alpha1'] == '-9.0'
    assert ki[('alprazolam', 'CHEMBL5143601')]['alpha2_minus_alpha1'] == '-0.2'
    assert len({p['integration_pair_id'] for p in pairs}) == 4  # context ID alone is shared by two drugs
    assert all(p['standard_units'] == 'nM' for p in pairs if p['ratio_eligible'])


@pytest.mark.parametrize('corruption', ['multiple_observations', 'compound_id', 'assay_id', 'species', 'duplicate_activity'])
def test_pair_identity_and_context_fail_closed(real, corruption):
    config, _ = real
    root = it.ROOT / config['pharmacology_run']
    pairs = it.read_csv(root / 'tables/alpha1_alpha2_pair_candidates.csv')
    rows = it.read_jsonl(root / 'processed/bzd_activity_classification.jsonl')
    if corruption == 'multiple_observations': pairs[0]['alpha1_activity_ids'] = '[667545,667547]'
    elif corruption == 'compound_id': pairs[0]['molecule_chembl_id'] = 'CHEMBL_FAKE'
    elif corruption == 'assay_id': pairs[0]['alpha1_assay_ids'] = '["CHEMBL_FAKE"]'
    elif corruption == 'duplicate_activity': rows.append(rows[0])
    else:
        next(r for r in rows if r['activity_id'] == 667545)['species'] = 'rat'
    with pytest.raises(ValueError): it.extract_pairs(pairs, rows, config['compound_ids'])


def structural_inputs(config):
    root = it.ROOT / config['structure_run']
    return (it.read_csv(root / 'data/processed/reassigned_interactions.csv'),
            it.read_csv(root / 'data/raw/lineage/docking_results_4drug.csv'),
            it.read_csv(root / 'results/tables/corrected_residue_mapping.csv'),
            it.read_csv(root / 'results/tables/candidate_residue_features.csv'),
            {sub: it.pdb_residues(it.ROOT / 'data/raw/structures/strict3' / (pdb + '.pdb')) for sub,pdb in [('alpha1','6HUP'),('alpha2','9CTJ')]})


def test_frequency_counts_poses_not_multiple_contacts(real):
    features = real[1][2]
    r = next(r for r in features if r['drug_id'] == 'alprazolam' and r['receptor_subtype'] == 'alpha1' and r['site_label'] == 'alpha_HIS102' and r['interaction_type'] == 'pi_stack')
    assert r['interaction_count'] == 2 and r['n_interacting_poses'] == 1
    assert r['frequency'] == 1/9 and r['type_T_frequency'] == 1/9
    assert r['centdist_n_observed'] == 2


def test_absent_pi_contacts_have_zero_frequency_but_missing_geometry(real):
    r = next(r for r in real[1][2] if r['drug_id'] == 'diazepam' and r['receptor_subtype'] == 'alpha2' and r['site_label'] == 'gamma2_PHE77' and r['interaction_type'] == 'pi_stack')
    assert r['frequency'] == 0 and r['centdist_mean'] is None and r['angle_mean'] is None
    assert r['feature_origin'] == 'requested_feature_audited_against_complete_raw_contacts'


@pytest.mark.parametrize('missing', ['pose_archive', 'residue'])
def test_missing_evidence_not_zero(real, missing):
    args = list(structural_inputs(real[0]))
    if missing == 'pose_archive':
        args[0] = [r for r in args[0] if not (r['drug_id'] == 'diazepam' and r['receptor'] == 'alpha2' and r['pose_id'] == '9')]
    else:
        args[4]['alpha2'].discard(('E', 77, 'PHE'))
    features, _, _ = it.fingerprint_rows(*args)
    r = next(r for r in features if r['drug_id'] == 'diazepam' and r['receptor_subtype'] == 'alpha2' and r['site_label'] == 'gamma2_PHE77' and r['interaction_type'] == 'pi_stack')
    assert r['frequency'] is None and r['data_status'] == 'missing_or_unresolved_evidence'
    _, matrix = it.heatmap_data(features)
    assert np.isnan(matrix).any()


def test_prior_frequency_tampering_detected(real):
    args = list(structural_inputs(real[0]))
    args[3][0]['frequency'] = '0.99'
    with pytest.raises(ValueError, match='frequency'): it.fingerprint_rows(*args)


def test_receptor_specific_actual_residue_numbers_and_all_features(real):
    features = real[1][2]
    assert len(features) == 48
    for r in features:
        if r['site_label'] == 'alpha_HIS102':
            assert r['actual_residue_number'] == (102 if r['receptor_subtype'] == 'alpha1' else 101)
        if r['site_label'] == 'gamma2_ASN60': assert r['interaction_type'] in {'hydrogen_bond', 'hydrophobic_interaction', 'halogen_bond'}


def test_integration_preserves_original_fields_and_flags(real):
    _, data = real
    activities, integrated = data[1], data[6]
    for a,r in zip(activities, integrated):
        assert all(r[k] == v for k,v in a.items())
        assert 'ifp__gamma2_TYR58__pi_stack__frequency' in r
        assert 'local_construct' in r['structure_context_limit']
        if 'beta2' in r['receptor_composition']: assert 'beta2_pharmacology_vs_beta3' in r['structure_context_limit']
    assert len(data[7]) == 48  # not 48 independent pharmacology observations


def test_geometry_keeps_raw_pt_contacts(real):
    geometry = real[1][4]
    p = [r for r in geometry if r['drug_id'] == 'alprazolam' and r['site_label'] == 'gamma2_PHE77']
    assert sorted(r['stacking_type'] for r in p) == ['P','P','P','T']
    assert len({r['source_row_id'] for r in geometry}) == len(geometry) == 23


def test_immutable_run_and_output_roundtrip(tmp_path, monkeypatch, real):
    original = it.ROOT
    # Symlink read-only input trees; output is isolated in tmp_path/runs.
    (tmp_path/'runs').mkdir(); (tmp_path/'runs/LATEST_RUN.txt').write_text('unchanged\n')
    for name in [real[0]['pharmacology_run'], real[0]['structure_run']]: (tmp_path/name).symlink_to(original/name, target_is_directory=True)
    (tmp_path/'data').symlink_to(original/'data', target_is_directory=True)
    (tmp_path/'src').symlink_to(original/'src', target_is_directory=True)
    (tmp_path/'tests').symlink_to(original/'tests', target_is_directory=True)
    shutil.copyfile(original/'pytest.ini', tmp_path/'pytest.ini')
    config = tmp_path/'config.json'; config.write_text(json.dumps(real[0]))
    pointer = tmp_path/'runs/LATEST_RUN.txt'
    monkeypatch.setattr(it, 'ROOT', tmp_path)
    monkeypatch.setattr(it.run_manager, 'ROOT', tmp_path)
    monkeypatch.setattr(it.run_manager, 'RUNS', tmp_path/'runs')
    monkeypatch.setattr(it, 'tracked_inventory', lambda: {'runs/LATEST_RUN.txt':it.sha(pointer.read_bytes())})
    # The actual renderer is exercised: images/PDFs must exist and be nonempty.
    out, qc = it.run(config)
    assert qc['status'] == 'PASS' and all(qc['checks'].values())
    assert pointer.read_text() == 'unchanged\n'
    assert len(it.read_csv(out/'matched_pharmacology_pairs.csv')) == 4
    assert len(it.read_csv(out/'diazepam_alprazolam_ifp.csv')) == 48
    m = json.loads((out/'run_manifest.json').read_text())
    assert all(it.sha((out/name).read_bytes()) == digest for name,digest in m['artifact_sha256'].items())
    assert len(list((out/'figures').glob('*.png'))) == 4
